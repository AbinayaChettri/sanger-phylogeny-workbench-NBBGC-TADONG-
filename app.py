import datetime
import gc
import math
import os
import platform
import re
import shutil
import subprocess
import tempfile
import threading
import zipfile
from collections import Counter, defaultdict

import Bio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from Bio import Align, Entrez, Phylo, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

APP_VERSION = "v2.1"
REPO_URL = "https://github.com/AbinayaChettri/sanger-phylogeny-workbench-NBBGC-TADONG-"
APP_URL = "https://nbbgc-tadong-phylogeny.streamlit.app"
AUTHOR_ORCID = "0009-0000-9037-9699"
FIGURE_CREDIT = f"Generated with Universal Sanger-to-Phylogeny Workbench ({APP_VERSION}) — Chettri (2026) | {REPO_URL}"

st.set_page_config(
    page_title=f"Sanger-to-Phylogeny Workbench ({APP_VERSION})",
    page_icon="🧬",
    layout="wide",
)


@st.cache_resource
def get_compute_queue():
    """Allows up to 2 simultaneous heavy alignment/IQ-TREE jobs to prevent OOM crashes on shared cloud servers."""
    return threading.BoundedSemaphore(2)


@st.cache_data(ttl=86400, show_spinner=False)
def get_tool_version(cmd_tuple, fallback="installed"):
    try:
        res = subprocess.run(list(cmd_tuple), capture_output=True, text=True, timeout=5)
        txt = (res.stdout + "\n" + res.stderr).strip().splitlines()
        return txt[0].strip()[:80] if txt else fallback
    except Exception:
        return fallback


def mott_trim(quals, q_cutoff=20):
    err_lim = 10 ** (-q_cutoff / 10.0)
    run_s = max_s = 0.0
    best_s = best_e = cur_s = 0
    for i, q in enumerate(quals):
        run_s += err_lim - 10 ** (-max(q, 1) / 10.0)
        if run_s <= 0:
            run_s = 0.0
            cur_s = i + 1
        elif run_s > max_s:
            max_s = run_s
            best_s, best_e = cur_s, i + 1
    return (best_s, best_e) if best_e > best_s else (0, len(quals))


def plot_qc(rec, s, e, out_png):
    quals = list(rec.letter_annotations.get("phred_quality", []))
    raw = rec.annotations.get("abif_raw", {})
    fig, axes = plt.subplots(2, 1, figsize=(11, 5.2), dpi=250, gridspec_kw={"height_ratios": [1.3, 1]})
    if all(k in raw for k in ("DATA9", "DATA10", "DATA11", "DATA12")):
        ploc = raw.get("PLOC2") or raw.get("PLOC1")
        s_idx = max(0, ploc[s] - 10) if ploc and len(ploc) > s else 0
        e_idx = min(len(raw["DATA9"]), ploc[min(s + 80, e - 1)] + 10) if ploc and len(ploc) > s else min(1000, len(raw["DATA9"]))
        xs = np.arange(s_idx, e_idx)
        for base, (k, col) in {"G": ("DATA9", "#111827"), "A": ("DATA10", "#059669"), "T": ("DATA11", "#DC2626"), "C": ("DATA12", "#2563EB")}.items():
            axes[0].plot(xs, raw[k][s_idx:e_idx], color=col, lw=0.85, label=base)
        axes[0].set_title(f"Chromatogram Trace Preview (Trimmed Start: Base {s+1})", fontweight="bold")
        axes[0].legend(loc="upper right", ncol=4, frameon=False)
    else:
        axes[0].text(0.5, 0.5, "FASTA input (no raw ABIF chromatogram channels)", ha="center", va="center", color="#4B5563")
        axes[0].set_axis_off()

    if quals:
        axes[1].fill_between(range(1, len(quals) + 1), quals, color="#0F766E", alpha=0.22)
        axes[1].plot(range(1, len(quals) + 1), quals, color="#0F766E", lw=1)
        axes[1].axhline(20, color="#DC2626", ls="--", lw=0.9, label="Q20 cutoff")
        axes[1].axvline(s + 1, color="#111827", ls=":", lw=1.1, label="Trim bounds")
        axes[1].axvline(e, color="#111827", ls=":", lw=1.1)
        axes[1].set_xlabel("Base position (bp)")
        axes[1].set_ylabel("Phred Q")
        axes[1].legend(frameon=False)
    else:
        axes[1].text(0.5, 0.5, "No Phred quality scores (FASTA input)", ha="center", va="center", color="#4B5563")
        axes[1].set_axis_off()
    fig.text(0.99, 0.008, FIGURE_CREDIT, ha="right", va="bottom", fontsize=6.5, color="#6B7280", style="italic")
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(out_png, dpi=250)
    plt.close(fig)


def merge_fr(seq_f, q_f, seq_r, q_r):
    """Aligns Forward and Reverse-Complement reads and returns consensus + overlap QC stats."""
    seq_r_rc = str(Seq(seq_r).reverse_complement())
    q_r_rc = list(reversed(q_r)) if q_r else []
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2.0, -1.0
    aligner.open_gap_score, aligner.extend_gap_score = -4.0, -0.5
    aln = aligner.align(seq_f, seq_r_rc)[0]
    a_str, b_str = str(aln[0]), str(aln[1])
    cons, i_f, i_r, mm, ov_bp, ov_match = [], 0, 0, 0, 0, 0
    for ca, cb in zip(a_str, b_str):
        qf = q_f[i_f] if (ca != "-" and i_f < len(q_f)) else -1
        qr = q_r_rc[i_r] if (cb != "-" and i_r < len(q_r_rc)) else -1
        if ca != "-":
            i_f += 1
        if cb != "-":
            i_r += 1
        if ca == "-":
            cons.append(cb)
        elif cb == "-":
            cons.append(ca)
        else:
            ov_bp += 1
            if ca == cb:
                ov_match += 1
                cons.append(ca)
            else:
                mm += 1
                cons.append(ca if qf >= qr else cb)
    ov_id = round(100.0 * ov_match / max(ov_bp, 1), 2) if ov_bp > 0 else 0.0
    return "".join(cons), mm, ov_bp, ov_id


def fetch_genbank_with_metadata(acc_list, outgroup_acc=None):
    """Fetches authoritative GenBank records (Organism, Marker, Country, Host) via Entrez."""
    refs_data, meta_rows = [], []
    outgroup_label = None
    with Entrez.efetch(db="nucleotide", id=",".join(acc_list), rettype="gb", retmode="text") as h:
        for gb in SeqIO.parse(h, "genbank"):
            acc = gb.id.split(".")[0]
            org = gb.annotations.get("organism", "Unknown")
            country, host, marker = "", "", ""
            for feat in gb.features:
                if feat.type == "source":
                    c_raw = feat.qualifiers.get("geo_loc_name", feat.qualifiers.get("country", [""]))[0]
                    country = c_raw.split(":")[0].strip()
                    host = feat.qualifiers.get("host", [""])[0].strip()
                elif not marker and feat.type in ("gene", "CDS", "rRNA", "misc_RNA"):
                    marker = feat.qualifiers.get("gene", feat.qualifiers.get("product", [""]))[0].strip()
            org_short = "_".join(org.split()[:2])
            parts = [acc, org_short]
            if country:
                parts.append(country)
            if host:
                parts.append(host.split()[0])
            clean_lbl = re.sub(r"[^A-Za-z0-9_.-]", "_", "_".join(parts))[:48]
            if outgroup_acc and acc.startswith(outgroup_acc.split(".")[0]):
                outgroup_label = clean_lbl
            refs_data.append((clean_lbl, str(gb.seq)))
            meta_rows.append({
                "Accession": acc,
                "Organism": org,
                "Marker": marker or "N/A",
                "Country": country or "N/A",
                "Host": host or "N/A",
                "Length_bp": len(gb.seq),
            })
    return refs_data, pd.DataFrame(meta_rows), outgroup_label


@st.cache_data(ttl=3600, max_entries=25, show_spinner=False)
def cached_blast_and_genbank(query_seq, max_refs, extra_accs_str, outgroup_acc, email):
    """Caches identical BLASTn + GenBank queries for 1 hour to reduce NCBI latency and server load."""
    Entrez.email = email.strip()
    handle = NCBIWWW.qblast("blastn", "nt", query_seq, hitlist_size=100)
    blast_rec = NCBIXML.read(handle)
    handle.close()

    all_hits, seen, blast_warnings = [], set(), []
    for aln in blast_rec.alignments:
        if aln.accession in seen:
            continue
        seen.add(aln.accession)
        hsp = aln.hsps[0]
        pct_id = round(100.0 * hsp.identities / max(hsp.align_length, 1), 2)
        q_cov = round(min(100.0, 100.0 * (abs(hsp.query_end - hsp.query_start) + 1) / max(len(query_seq), 1)), 1)
        title = aln.title.split("|")[-1].strip()
        words = title.split()
        species = " ".join(words[:2]) if len(words) >= 2 else title
        all_hits.append({"Accession": aln.accession, "Species": species, "Title": title[:90], "Pct_ID": pct_id, "Cov_Pct": q_cov, "E_Value": hsp.expect})

    blast_df = pd.DataFrame(all_hits)

    top_sp = {r["Species"] for r in all_hits if r["Pct_ID"] >= 99.0 and r["Cov_Pct"] >= 80.0}
    if len(top_sp) > 1:
        blast_warnings.append(f"BLAST multi-species tie: multiple taxa matched >=99% ID ({', '.join(sorted(top_sp)[:4])}). Inspect tree & metadata carefully.")
    if all_hits and all_hits[0]["Pct_ID"] < 97.0:
        blast_warnings.append(f"Low top BLAST identity ({all_hits[0]['Pct_ID']}%). Species-level ID is uncertain; verify marker and chromatogram.")

    by_sp = defaultdict(list)
    for row in all_hits:
        if row["Cov_Pct"] >= 55.0:
            by_sp[row["Species"]].append(row)

    picked, picked_accs = [], set()
    sp_order = list(by_sp.keys())
    if sp_order:
        for r in by_sp[sp_order[0]][:4]:
            picked.append(r)
            picked_accs.add(r["Accession"])
        for sp in sp_order[1:]:
            for r in by_sp[sp][:3]:
                if len(picked) < max_refs:
                    picked.append(r)
                    picked_accs.add(r["Accession"])
    for r in all_hits:
        if len(picked) >= max_refs:
            break
        if r["Accession"] not in picked_accs and r["Cov_Pct"] >= 55.0:
            picked.append(r)
            picked_accs.add(r["Accession"])

    acc_list = [r["Accession"] for r in picked]
    if extra_accs_str and extra_accs_str.strip():
        for a in re.split(r"[\s,;]+", extra_accs_str.strip()):
            if a and a not in acc_list:
                acc_list.append(a)

    og_target = outgroup_acc.strip() if outgroup_acc and outgroup_acc.strip() else None
    if og_target and og_target not in acc_list:
        acc_list.append(og_target)

    refs_data, meta_df, outgroup_label = fetch_genbank_with_metadata(acc_list, og_target)
    if og_target and not outgroup_label:
        blast_warnings.append(f"Requested outgroup '{og_target}' could not be retrieved from GenBank; falling back to midpoint rooting.")

    for _, row in meta_df.iterrows():
        if row["Length_bp"] < 0.60 * len(query_seq):
            blast_warnings.append(f"Short reference/outgroup detected: {row['Accession']} ({row['Organism']}, {row['Length_bp']} bp vs. query {len(query_seq)} bp) — short sequences can distort branch lengths.")

    return refs_data, blast_df, meta_df, outgroup_label, blast_warnings


def smart_blast_v2(query_seq, max_refs, extra_accs_str, outgroup_acc, email, work_dir, warnings_list):
    refs_data, blast_df, meta_df, outgroup_label, blast_warnings = cached_blast_and_genbank(
        query_seq, max_refs, extra_accs_str or "", outgroup_acc or "", email.strip()
    )
    warnings_list.extend(blast_warnings)
    blast_df.to_csv(f"{work_dir}/02_blast_top_hits.csv", index=False)
    meta_df.to_csv(f"{work_dir}/02_reference_metadata_country_host.csv", index=False)
    refs = [SeqRecord(Seq(seq_s), id=lbl, description="") for lbl, seq_s in refs_data]
    return refs, blast_df.head(20), meta_df, outgroup_label


def compute_popgen_stats(trimmed_recs, work_dir, outgroup_label=None):
    """Computes PopGen stats on INGROUP sequences only (excluding biological outgroup) and reports excluded gapped sites."""
    ingroup = [r for r in trimmed_recs if not (outgroup_label and r.id == outgroup_label)]
    if len(ingroup) < 2:
        ingroup = trimmed_recs

    seqs = [str(r.seq).upper() for r in ingroup]
    n = len(seqs)
    L = len(seqs[0])
    arr = np.array([list(s) for s in seqs])

    clean_cols = [c for c in range(L) if all(b in "ACGT" for b in arr[:, c])]
    excluded_cols = L - len(clean_cols)
    sub = arr[:, clean_cols] if clean_cols else arr
    L_clean = sub.shape[1]

    hap_strs = ["".join(sub[i]) for i in range(n)]
    hap_counts = Counter(hap_strs)
    h = len(hap_counts)
    hd = (n / (n - 1.0)) * (1.0 - sum((c / n) ** 2 for c in hap_counts.values())) if n > 1 else 0.0

    S = sum(1 for c in range(L_clean) if len(set(sub[:, c])) > 1)
    pi_sites = sum(1 for c in range(L_clean) if sum(1 for v in Counter(sub[:, c]).values() if v >= 2) >= 2)
    pair_diffs = [sum(sub[i, c] != sub[j, c] for c in range(L_clean)) for i in range(n) for j in range(i + 1, n)]
    k = float(np.mean(pair_diffs)) if pair_diffs else 0.0
    pi_val = k / max(L_clean, 1)

    tajima_d = float("nan")
    if n >= 4 and S > 0:
        a1 = sum(1.0 / i for i in range(1, n))
        a2 = sum(1.0 / (i * i) for i in range(1, n))
        b1 = (n + 1.0) / (3.0 * (n - 1.0))
        b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
        c1 = b1 - 1.0 / a1
        c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1)
        e1 = c1 / a1
        e2 = c2 / (a1 * a1 + a2)
        var_d = e1 * S + e2 * S * (S - 1.0)
        if var_d > 0:
            tajima_d = (k - S / a1) / math.sqrt(var_d)

    pop_df = pd.DataFrame([{
        "Ingroup_Taxa_n": int(n),
        "Ungapped_ACGT_Sites_bp": int(L_clean),
        "Excluded_Gapped_or_Ambig_bp": int(excluded_cols),
        "Haplotypes_h": int(h),
        "Haplotype_Diversity_Hd": round(hd, 4),
        "Nucleotide_Diversity_pi": round(pi_val, 5),
        "Segregating_Sites_S": int(S),
        "Parsimony_Informative": int(pi_sites),
        "Tajimas_D": round(tajima_d, 4) if not math.isnan(tajima_d) else "N/A",
    }])
    pop_df.to_csv(f"{work_dir}/05_population_genetics_summary.csv", index=False)
    return pop_df


def plot_snp_map(trimmed_recs, out_png):
    arr = np.array([list(str(r.seq).upper()) for r in trimmed_recs])
    var_cols = [c for c in range(arr.shape[1]) if len(set(arr[:, c]) - {"-", "N"}) > 1]
    fig, ax = plt.subplots(figsize=(11, max(4.2, len(trimmed_recs) * 0.35 + 1.4)), dpi=220)
    if not var_cols:
        ax.text(0.5, 0.5, "No variable SNP sites detected across trimmed alignment.", ha="center", va="center")
        ax.set_axis_off()
    else:
        show_cols = var_cols[:60]
        cmap = {"A": "#059669", "T": "#DC2626", "C": "#2563EB", "G": "#F59E0B", "-": "#E5E7EB", "N": "#9CA3AF"}
        for r_idx, rec in enumerate(trimmed_recs):
            for x_idx, col in enumerate(show_cols):
                b = arr[r_idx, col]
                ax.add_patch(plt.Rectangle((x_idx, len(trimmed_recs) - 1 - r_idx), 1, 1, facecolor=cmap.get(b, "#D1D5DB"), edgecolor="white", lw=0.4))
                ax.text(x_idx + 0.5, len(trimmed_recs) - 0.5 - r_idx, b, ha="center", va="center", fontsize=6.5, color="white" if b in "ATCG" else "#111827", fontweight="bold")
        ax.set_xlim(0, len(show_cols))
        ax.set_ylim(0, len(trimmed_recs))
        ax.set_yticks([i + 0.5 for i in range(len(trimmed_recs))])
        ax.set_yticklabels([r.id[:32] for r in reversed(trimmed_recs)], fontsize=7.5)
        ax.set_xticks([i + 0.5 for i in range(len(show_cols))])
        ax.set_xticklabels([str(c + 1) for c in show_cols], rotation=90, fontsize=6.5)
        ax.set_xlabel("Alignment Position (bp) — All Variable Sites (including gapped/outgroup columns)")
        ax.set_title(f"Alignment SNP & Mismatch Map ({len(var_cols)} variable sites; showing first {len(show_cols)})", fontweight="bold", fontsize=10)
    fig.text(0.99, 0.008, FIGURE_CREDIT, ha="right", va="bottom", fontsize=6.5, color="#6B7280", style="italic")
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(out_png, dpi=220)
    plt.close(fig)
    return len(var_cols)


# ==================== STREAMLIT UI ====================
st.title(f"🧬 Universal Sanger-to-Phylogeny Workbench ({APP_VERSION})")
st.caption(f"Developed by **Abinaya Chettri** (M.Sc. Zoology | ORCID: {AUTHOR_ORCID}) — Automated Sanger QC, Validated F+R Consensus, NCBI BLASTn + Country/Host Metadata, MAFFT, SNP Map, Ingroup Population Genetics & IQ-TREE 2 Phylogenetics.")

with st.expander("📖 Quick-Start Rules & Outgroup Guide (Read First Before Running)", expanded=False):
    st.markdown("""
**1. `.ab1` vs. `.fasta` Rule:**
* For your own PCR samples from a sequencing company, **always upload Forward `.ab1` (Box 1) + Reverse `.ab1` (Box 2)**—they hold raw trace peaks and Phred quality scores.
* Use `.fasta` only for already-cleaned consensus sequences or published GenBank records.

**2. Two-Run Workflow (Discovery → Publication Figure):**
* **Run 1 (Discovery):** Keep Phred Cutoff at `20` and leave both Optional boxes **blank**. Let BLASTn identify your species and exact sequence length (e.g., `658 bp` of `cox1`).
* **Run 2 (Paper Figure):** Paste 2–4 verified GenBank accessions from published papers + 1 **Biological Outgroup** accession.

**3. Outgroup Golden Rules (Prevents Squished Trees):**
* **Same gene & similar length:** Match `cox1` with `cox1` (or `ITS-2` with `ITS-2`) of roughly the same length (≥70–100% of your query length).
* **Close cousin, not a distant stranger:** Choose a sister species in the **same genus** (e.g., *Toxocara cati* for *Toxocara canis*; *Ancylostoma caninum* for *Ancylostoma ceylanicum*).
* **Strong clade support threshold:** SH-aLRT ≥ 80% and UFBoot ≥ 95%.

**4. Population Genetics & Publication Note:**
* v2.1 automatically **excludes the Outgroup** when computing `h`, `Hd`, `π`, `S`, and Tajima's `D`, and uses **strict ungapped A/C/G/T sites** (standard DnaSP/Tajima convention). For formal demographic inference (Tajima's `D`), upload your own single-population sample set via Box 3 rather than mixed global BLAST hits.
* Upload `01_clean_query.fasta` to **NCBI Submission Portal — GenBank** for official accession numbers, and open `06_tree_for_MEGA_FigTree_iTOL.nwk` in **iTOL** (`itol.embl.de`) or **FigTree** for custom layouts.
""")

col1, col2, col3 = st.columns(3)
with col1:
    primary_up = st.file_uploader("1. Upload Query (.ab1 or .fasta) *")
with col2:
    reverse_up = st.file_uploader("2. Optional Reverse Read (.ab1 for F+R consensus)")
with col3:
    custom_refs_up = st.file_uploader("3. Optional Custom References (.fasta — skips BLAST)")

c1, c2 = st.columns(2)
with c1:
    q_cutoff = st.slider("Phred Trim Cutoff (Mott Q)", 10, 40, 20, 1)
with c2:
    max_refs = st.slider("Total Diverse BLAST References", 5, 30, 14, 1)

c3, c4, c5 = st.columns(3)
with c3:
    extra_accs = st.text_input("Optional Published Reference Accessions (comma-separated)", placeholder="e.g., LC533328, AJ407962")
with c4:
    outgroup_acc = st.text_input("Optional Biological Outgroup Accession", placeholder="e.g., AJ407962 (blank = midpoint)")
with c5:
    email = st.text_input("NCBI Email (required for live BLAST) *", value=os.environ.get("NCBI_EMAIL", ""), placeholder="your.email@institution.edu")

if st.button(f"🚀 Run Complete {APP_VERSION} Pipeline", type="primary", use_container_width=True):
    if primary_up is None:
        st.error("Please upload a `.ab1` chromatogram or `.fasta` file in Box 1 first.")
        st.stop()
    if custom_refs_up is None and ("@" not in email.strip()):
        st.error("Please enter a valid email address in the NCBI Email box (NCBI Entrez requires it for live BLAST/GenBank lookups).")
        st.stop()

    iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")
    if shutil.which("mafft") is None or iqtree_bin is None:
        st.error("Missing system binaries (`mafft` or `iqtree`). Ensure `packages.txt` includes `mafft` and `iqtree`.")
        st.stop()

    work = tempfile.mkdtemp(prefix="phylo_v21_")
    notes, warnings_list = [], []
    prog = st.progress(5, text="Reading & quality-trimming input...")

    try:
        def save_upload(up_file, subdir):
            p = os.path.join(subdir, up_file.name)
            with open(p, "wb") as f:
                f.write(up_file.getbuffer())
            return p

        def read_seq(path):
            fmt = "abi" if path.lower().endswith((".ab1", ".abi")) else "fasta"
            rec = next(SeqIO.parse(path, fmt))
            raw_s = str(rec.seq).upper()
            quals = list(rec.letter_annotations.get("phred_quality", []))
            if quals:
                s, e = mott_trim(quals, int(q_cutoff))
                trim_s, trim_q = raw_s[s:e], quals[s:e]
                mean_q = float(np.mean(trim_q)) if trim_q else 0.0
                ambig = sum(1 for b in trim_s if b not in "ACGT")
                notes.append(f"{os.path.basename(path)}: {len(raw_s)} bp -> {len(trim_s)} bp (mean Q={mean_q:.1f}, ambiguous={ambig})")
                if len(trim_s) < 150:
                    warnings_list.append(f"{os.path.basename(path)} trimmed to only {len(trim_s)} bp (<150 bp); check chromatogram quality.")
                if mean_q < 25.0:
                    warnings_list.append(f"{os.path.basename(path)} has low mean quality (Q={mean_q:.1f}); inspect peaks manually.")
            else:
                s, e = 0, len(raw_s)
                trim_s, trim_q = re.sub(r"[^ACGTNRYMKSWBDHV]", "", raw_s), []
                notes.append(f"{os.path.basename(path)} (FASTA): {len(trim_s)} bp")
            return rec, trim_s, trim_q, s, e

        p_path = save_upload(primary_up, work)
        rec_f, seq_f, q_f, s_f, e_f = read_seq(p_path)
        qc_png = f"{work}/01_chromatogram_qc.png"
        plot_qc(rec_f, s_f, e_f, qc_png)

        final_seq = seq_f
        if reverse_up is not None:
            r_path = save_upload(reverse_up, work)
            _, seq_r, q_r, _, _ = read_seq(r_path)
            final_seq, mm, ov_bp, ov_id = merge_fr(seq_f, q_f, seq_r, q_r)
            notes.append(f"F+R Consensus: {len(final_seq)} bp | Overlap={ov_bp} bp ({ov_id}% identity, {mm} conflicts resolved by higher Phred Q)")
            if ov_bp < 50 or ov_id < 90.0:
                warnings_list.append(f"Low Forward/Reverse consensus confidence (overlap={ov_bp} bp, identity={ov_id}%). Verify both reads belong to the same PCR amplicon before trusting results!")

        q_id = "QUERY_" + re.sub(r"[^A-Za-z0-9_.-]", "_", rec_f.id or "Sample")[:28]
        query_rec = SeqRecord(Seq(final_seq), id=q_id, description="")
        SeqIO.write([query_rec], f"{work}/01_clean_query.fasta", "fasta")

        meta_df = pd.DataFrame()
        if custom_refs_up is not None:
            prog.progress(30, text="Loading custom reference FASTA...")
            c_path = save_upload(custom_refs_up, work)
            refs = list(SeqIO.parse(c_path, "fasta"))
            blast_df = pd.DataFrame([{"Accession": r.id, "Length_bp": len(r.seq)} for r in refs])
            outgroup_label = outgroup_acc.strip() if outgroup_acc and outgroup_acc.strip() else None
        else:
            prog.progress(20, text="Running live NCBI BLASTn + fetching GenBank metadata (~30-60s)...")
            refs, blast_df, meta_df, outgroup_label = smart_blast_v2(final_seq, int(max_refs), extra_accs, outgroup_acc, email, work, warnings_list)

        if len(refs) < 3:
            st.error("Need at least 3 reference sequences to build a phylogenetic tree.")
            st.stop()

        prog.progress(50, text="Waiting for server compute slot (keeps multi-user runs smooth)...")
        with get_compute_queue():
            prog.progress(55, text="Running MAFFT alignment & gap trimming...")
            raw_fa, aln_fa, trim_fa = f"{work}/03_unaligned.fasta", f"{work}/03_mafft_aligned.fasta", f"{work}/04_trimmed_aligned.fasta"
            seen_ids, uniq = {}, []
            for r in [query_rec] + refs:
                bid = re.sub(r"[^A-Za-z0-9_.-]", "_", r.id)[:48]
                seen_ids[bid] = seen_ids.get(bid, 0) + 1
                r.id = bid if seen_ids[bid] == 1 else f"{bid}_{seen_ids[bid]}"
                r.description = ""
                uniq.append(r)
            SeqIO.write(uniq, raw_fa, "fasta")
            with open(aln_fa, "w") as out_f:
                subprocess.run(["mafft", "--auto", "--adjustdirection", "--thread", "2", raw_fa], stdout=out_f, stderr=subprocess.DEVNULL, check=True)
            aligned = []
            for r in SeqIO.parse(aln_fa, "fasta"):
                if r.id.startswith("_R_"):
                    r.id = r.id[3:]
                r.description = ""
                aligned.append(r)
            arr = np.array([list(str(r.seq).upper()) for r in aligned])
            keep = [c for c in range(arr.shape[1]) if np.mean(arr[:, c] == "-") <= 0.5]
            trimmed = [SeqRecord(Seq("".join(arr[i, keep])), id=aligned[i].id, description="") for i in range(len(aligned))]
            SeqIO.write(trimmed, trim_fa, "fasta")

            snp_png = f"{work}/04_alignment_snp_map.png"
            total_snp_sites = plot_snp_map(trimmed, snp_png)

            prog.progress(68, text="Computing K2P distances & Ingroup Population Genetics...")
            pop_df = compute_popgen_stats(trimmed, work, outgroup_label=outgroup_label)
            names, seqs = [r.id for r in trimmed], [str(r.seq) for r in trimmed]
            n = len(seqs)
            p_mat, k_mat = np.zeros((n, n)), np.zeros((n, n))
            trans = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
            for i in range(n):
                for j in range(i + 1, n):
                    ti = tv = v = 0
                    for a, b in zip(seqs[i], seqs[j]):
                        if a in "ACGT" and b in "ACGT":
                            v += 1
                            if a != b:
                                if (a, b) in trans:
                                    ti += 1
                                else:
                                    tv += 1
                    P, Q = ti / max(v, 1), tv / max(v, 1)
                    t1, t2 = 1 - 2 * P - Q, 1 - 2 * Q
                    k2p = (-0.5 * math.log(t1) - 0.25 * math.log(t2)) if (v > 0 and t1 > 0 and t2 > 0) else float("nan")
                    p_mat[i, j] = p_mat[j, i] = abs(round((ti + tv) / max(v, 1), 5))
                    k_mat[i, j] = k_mat[j, i] = abs(round(k2p, 5)) if not math.isnan(k2p) else float("nan")
            pd.DataFrame(p_mat, index=names, columns=names).to_csv(f"{work}/05_p_distance_matrix.csv")
            pd.DataFrame(k_mat, index=names, columns=names).to_csv(f"{work}/05_k2p_distance_matrix.csv")

            prog.progress(80, text="Running IQ-TREE 2 (ModelFinder + 1,000x UFBoot + SH-aLRT)...")
            prefix = f"{work}/06_iqtree"
            iq_cmd = [iqtree_bin, "-s", trim_fa, "-m", "MFP", "-B", "1000", "--alrt", "1000", "-bnni", "-nt", "2", "--prefix", prefix, "--redo"]
            if outgroup_label and outgroup_label in names:
                iq_cmd.extend(["-o", outgroup_label])
            subprocess.run(iq_cmd, check=True, capture_output=True)

        shutil.copy(f"{prefix}.treefile", f"{work}/06_tree_for_MEGA_FigTree_iTOL.nwk")
        with open(f"{prefix}.iqtree") as rep_f:
            rep_txt = rep_f.read()
        m = re.search(r"Best-fit model according to BIC:\s*(\S+)", rep_txt)
        best_model = m.group(1) if m else "ModelFinder"

        prog.progress(92, text="Rendering 300 DPI tree & packaging ZIP...")
        tree = Phylo.read(f"{prefix}.treefile", "newick")
        root_note = "Midpoint-rooted (no biological outgroup specified)"
        if outgroup_label and outgroup_label in names:
            try:
                tree.root_with_outgroup(outgroup_label)
                root_note = f"Rooted on biological outgroup: {outgroup_label}"
            except Exception:
                tree.root_at_midpoint()
        else:
            tree.root_at_midpoint()
        tree.ladderize()

        fig_h = max(6.2, len(tree.get_terminals()) * 0.45 + 1.8)
        fig, ax = plt.subplots(figsize=(12.5, fig_h), dpi=300)
        Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: f"  {c.name}" if c.is_terminal() else "", show_confidence=False)
        max_x = max(tree.depths().values()) or 0.01
        ax.set_xlim(-max_x * 0.03, max_x * 1.42)
        for t in ax.texts:
            if q_id in t.get_text():
                t.set_color("#0F766E")
                t.set_fontweight("bold")
                t.set_fontsize(9.5)
            else:
                t.set_color("#1F2937")
                t.set_fontsize(8.5)

        depths = tree.depths()
        y_pos = {term: idx for idx, term in enumerate(tree.get_terminals(), start=1)}
        for clade in tree.get_nonterminals(order="postorder"):
            ys = [y_pos[c] for c in clade.clades if c in y_pos]
            if ys:
                y_pos[clade] = float(np.mean(ys))
                if clade.name and "/" in str(clade.name):
                    ax.text(depths.get(clade, 0.0) - max_x * 0.008, y_pos[clade] - 0.14, str(clade.name), fontsize=7.5, color="#374151", ha="right", va="bottom")

        ax.set_title(f"Maximum Likelihood Phylogeny (IQ-TREE 2 | Model: {best_model} | SH-aLRT / UFBoot 1000x)\n{root_note}", fontweight="bold", fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_xlabel("Substitutions per site")
        fig.text(0.99, 0.008, FIGURE_CREDIT, ha="right", va="bottom", fontsize=6.5, color="#6B7280", style="italic")
        plt.tight_layout(rect=[0, 0.025, 1, 1])
        tree_png, tree_pdf = f"{work}/06_phylogenetic_tree_300dpi.png", f"{work}/06_phylogenetic_tree_vector.pdf"
        fig.savefig(tree_png, dpi=300)
        fig.savefig(tree_pdf)
        plt.close(fig)

        mafft_ver = get_tool_version(("mafft", "--version"))
        iqtree_ver = get_tool_version((iqtree_bin, "--version"))
        acc_used = ", ".join(meta_df["Accession"].tolist()) if not meta_df.empty else "Custom FASTA references"
        pr = pop_df.iloc[0]
        warn_block = "\n  - ".join(warnings_list) if warnings_list else "None"

        methods_citation = (
            f"Sanger sequence quality control (Mott trimming at Q={int(q_cutoff)}), NCBI BLASTn reference curation, "
            f"multiple sequence alignment (MAFFT --auto across {len(trimmed)} taxa; {len(keep)} bp after <=50% gap trimming), "
            f"ingroup population genetics, and Maximum-Likelihood phylogenetic inference (IQ-TREE 2; best-fit BIC model {best_model}; "
            f"1,000 SH-aLRT and 1,000 ultrafast bootstrap replicates; {root_note}) were performed using Universal Sanger-to-Phylogeny Workbench "
            f"({APP_VERSION}; Chettri, 2026; {REPO_URL}).\n\n"
            f"Reference:\n"
            f"Chettri, A. (2026). Universal Sanger-to-Phylogeny Workbench ({APP_VERSION}) [Computer software]. "
            f"ORCID: https://orcid.org/{AUTHOR_ORCID}. Available at: {REPO_URL} (Web app: {APP_URL})."
        )

        summary = (
            f"Universal Sanger-to-Phylogeny Workbench ({APP_VERSION}) — Developed by Abinaya Chettri (ORCID: {AUTHOR_ORCID})\n"
            f"Run Date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n"
            f"Software Versions: Python {platform.python_version()} | Biopython {Bio.__version__} | {mafft_ver} | {iqtree_ver}\n"
            f"Parameters: Mott Q={int(q_cutoff)} | Max BLAST refs={int(max_refs)} | Gap trim <=50%\n\n"
            f"HOW TO CITE (Ready-to-Paste Methods & Reference):\n"
            f"{methods_citation}\n\n"
            f"QC & Consensus:\n  - " + "\n  - ".join(notes) + "\n\n"
            f"Alignment & SNP Consistency:\n"
            f"  - Total taxa: {len(trimmed)} | Alignment: {arr.shape[1]} bp raw -> {len(keep)} bp trimmed\n"
            f"  - Visual SNP Map variable sites (all taxa, incl. gapped/outgroup cols): {total_snp_sites}\n"
            f"  - Ingroup PopGen strict ungapped ACGT sites: {int(pr['Ungapped_ACGT_Sites_bp'])} bp (excluded {int(pr['Excluded_Gapped_or_Ambig_bp'])} gapped/ambiguous cols)\n\n"
            f"Ingroup Population Genetics (Outgroup excluded):\n"
            f"  - n={int(pr['Ingroup_Taxa_n'])} | h={int(pr['Haplotypes_h'])} | Hd={pr['Haplotype_Diversity_Hd']} | pi={pr['Nucleotide_Diversity_pi']} | S={int(pr['Segregating_Sites_S'])} | Tajima's D={pr['Tajimas_D']}\n"
            f"  - Note: Interpret Tajima's D cautiously if ingroup consists of mixed global GenBank BLAST hits rather than a single biological population.\n\n"
            f"Phylogeny:\n"
            f"  - Best-fit Model (BIC): {best_model} | Support: 1,000x SH-aLRT / 1,000x UFBoot (-bnni)\n"
            f"  - Rooting: {root_note}\n"
            f"  - GenBank Accessions Used: {acc_used}\n\n"
            f"Validation Warnings:\n  - {warn_block}\n"
        )
        with open(f"{work}/07_methods_and_reproducibility_report.txt", "w") as f:
            f.write(summary)

        zip_name = f"phylogeny_v21_results_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
        zip_path = os.path.join(work, zip_name)
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for fn in sorted(os.listdir(work)):
                if fn.endswith(".zip"):
                    continue
                zf.write(os.path.join(work, fn), arcname=fn)

        with open(zip_path, "rb") as f:
            zip_bytes = f.read()
        with open(tree_png, "rb") as f:
            tree_bytes = f.read()
        with open(snp_png, "rb") as f:
            snp_bytes = f.read()
        with open(qc_png, "rb") as f:
            qc_bytes = f.read()

        st.session_state["last_run"] = {
            "warnings": warnings_list,
            "zip_name": zip_name,
            "zip_bytes": zip_bytes,
            "methods_citation": methods_citation,
            "tree_bytes": tree_bytes,
            "snp_bytes": snp_bytes,
            "qc_bytes": qc_bytes,
            "pop_df": pop_df,
            "blast_df": blast_df,
            "summary": summary,
        }
        prog.progress(100, text="Done!")
    finally:
        plt.close("all")
        shutil.rmtree(work, ignore_errors=True)
        gc.collect()

if "last_run" in st.session_state:
    res = st.session_state["last_run"]
    for w in res["warnings"]:
        st.warning(f"⚠️ {w}")
    st.success("Pipeline complete! Download your full reproducibility archive below (includes original uploaded files + logs).")
    st.download_button(
        f"⬇️ Download Complete {APP_VERSION} Results (.zip)",
        res["zip_bytes"],
        file_name=res["zip_name"],
        mime="application/zip",
        use_container_width=True,
    )

    st.info("📚 **Using these results in a thesis or paper?** Copy-paste the ready-made Methods & Citation text below:")
    st.code(res["methods_citation"], language="text")

    st.subheader("ML Phylogenetic Tree (300 DPI, Country + Host labels)")
    st.image(res["tree_bytes"], use_container_width=True)
    st.subheader("Colored Alignment SNP & Mismatch Map")
    st.image(res["snp_bytes"], use_container_width=True)
    st.subheader("Chromatogram & Phred QC")
    st.image(res["qc_bytes"], use_container_width=True)
    st.subheader("Ingroup Population Genetics & Haplotype Diversity (Outgroup Excluded)")
    st.caption("Computed on ingroup sequences using strict ungapped A/C/G/T sites (standard DnaSP/Tajima convention). Use Box 3 with a single biological population for formal demographic inference.")
    st.dataframe(res["pop_df"], use_container_width=True)
    st.subheader("Top BLAST Hits")
    st.dataframe(res["blast_df"], use_container_width=True)
    st.subheader("QC, Model & Reproducibility Report")
    st.code(res["summary"], language="text")
