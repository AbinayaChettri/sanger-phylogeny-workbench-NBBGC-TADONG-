# Universal Sanger-to-Phylogeny Workbench (v2.0)

**Developed by Abinaya Chettri (M.Sc. Zoology)**

An automated, zero-code browser workbench that turns raw Sanger `.ab1` chromatograms (Forward + Reverse) or `.fasta` sequences into publication-ready Maximum Likelihood phylogenies, colored SNP maps, genetic distance matrices, and population genetics summaries.

## What v2.0 Automates
1. **Chromatogram QC & Consensus:** Mott Q20 quality trimming, fluorescence trace + Phred Q plot, and Forward + Reverse (F+R) consensus building with Phred-weighted conflict resolution.
2. **Live Species ID & Metadata:** Automated NCBI BLASTn (`nt`) with automatic GenBank **Country + Host** metadata extraction directly onto tree tip labels.
3. **Alignment & SNP Visualization:** `MAFFT --auto --adjustdirection`, `<=50%` gap trimming, and a colored variable-site SNP & Mismatch Map.
4. **Genetic Distances & Population Genetics:** Pairwise `p-distance` and Kimura 2-Parameter (`K2P`) CSV matrices plus `h`, `Hd`, `pi`, `S`, parsimony-informative sites, and Tajima's `D`.
5. **Maximum Likelihood Phylogeny:** `IQ-TREE 2` with automatic `ModelFinder` (`-m MFP`), 1,000x Ultrafast Bootstrap (`-B 1000 -bnni`), and 1,000x SH-aLRT (`--alrt 1000`), exporting 300 DPI PNG, vector PDF, `.nwk` (for FigTree/iTOL), and a draft Methods report in one `.zip`.

## Quick-Start Rules (Run 1 vs. Run 2 & Outgroup Guide)
* **Always use `.ab1` (F + R) for your own PCR samples** from a sequencing company (`.fasta` is only for already-cleaned or GenBank sequences).
* **Keep Mott Phred Cutoff at `20`** and check `01_chromatogram_qc.png` for crisp single peaks.
* **Run 1 (Discovery):** Leave both Optional boxes blank so BLASTn identifies your species and exact gene length (e.g., `658 bp` of `cox1`).
* **Run 2 (Paper Figure):** Add 2–4 verified GenBank accessions from peer-reviewed papers + 1 **Biological Outgroup** following two golden rules:
  1. **Same gene & similar length** (at least 70–100% as long as your query—never mix a short 300 bp fragment into a 658 bp alignment).
  2. **Close sister species in the same genus** (e.g., *Toxocara cati* for *Toxocara canis*, or *Ancylostoma caninum* for *Ancylostoma ceylanicum*).
* **Publication Steps:** Upload `01_clean_query.fasta` (from the `.zip`) to the [NCBI Submission Portal — GenBank](https://submit.ncbi.nlm.nih.gov/about/genbank/) for official accession numbers, and open `06_tree_for_MEGA_FigTree_iTOL.nwk` in [iTOL](https://itol.embl.de/) or [FigTree](https://github.com/rambaut/figtree/releases) for custom layouts.

## Free 1-Click Deployment on Streamlit Community Cloud
1. Create a public GitHub repository named `sanger-phylogeny-workbench`.
2. Upload all 6 files (`app.py`, `requirements.txt`, `packages.txt`, `README.md`, `LICENSE`, `CITATION.cff`).
3. Sign in to Streamlit Community Cloud with GitHub, select `sanger-phylogeny-workbench` -> `app.py`, and click **Deploy**.

## License & Citation
Released under the MIT License — Copyright (c) 2026 Abinaya Chettri. Click **"Cite this repository"** on GitHub for the formal citation metadata (`CITATION.cff`).
