# Universal Sanger-to-Phylogeny Workbench (v2.1)

**Developed by Abinaya Chettri** (M.Sc. Zoology | ORCID: [0009-0000-9037-9699](https://www.google.com/search?q=https%3A%2F%2Forcid.org%2F0009-0000-9037-9699))  
**Live Web App (Free, Zero Install):** [https://nbbgc-tadong-phylogeny.streamlit.app](https://nbbgc-tadong-phylogeny.streamlit.app)

An automated, zero-code browser workbench that turns raw Sanger `.ab1` chromatograms (Forward + Reverse) or clean `.fasta` sequences into publication-ready Maximum-Likelihood phylogenies, colored SNP maps, genetic distance matrices, and ingroup population genetics summaries—all exported in one reproducible `.zip` archive.

---

## What v2.1 Automates

* **Chromatogram QC & Validated F+R Consensus:** Mott Q20 Phred quality trimming, fluorescence trace + Phred quality plot (`01_chromatogram_qc.png`), and Phred-weighted Forward + Reverse consensus building with overlap/identity QC warnings (`01_clean_query.fasta`).
* **Live Species ID & GenBank Metadata:** Automated NCBI BLASTn (`nt`) search with automatic GenBank Country + Host metadata extraction directly onto tree tip labels (`02_blast_top_hits.csv`, `02_reference_metadata_country_host.csv`).
* **Alignment & SNP Visualization:** `MAFFT --auto --adjustdirection`, `<=50%` gap trimming, and a colored variable-site SNP & Mismatch Map (`04_alignment_snp_map.png`).
* **Genetic Distances & Ingroup Population Genetics:** Pairwise p-distance and Kimura 2-Parameter (K2P) CSV matrices plus ingroup-only (outgroup-excluded) haplotype count (`h`), haplotype diversity (`Hd`), nucleotide diversity (`π`), segregating sites (`S`), parsimony-informative sites, and Tajima's `D` on strict ungapped A/C/G/T sites (`05_population_genetics_summary.csv`).
* **Maximum-Likelihood Phylogeny:** IQ-TREE 2 with automatic ModelFinder (`-m MFP`), 1,000× Ultrafast Bootstrap (`-B 1000 -bnni`), and 1,000× SH-aLRT (`--alrt 1000`), exporting a 300 DPI PNG, vector PDF, Newick tree (`06_tree_for_MEGA_FigTree_iTOL.nwk` for FigTree/iTOL/MEGA), and a ready-to-paste Methods & Reproducibility report (`07_methods_and_reproducibility_report.txt`).

---

## Multi-Phylum Validation

Tested end-to-end across three independent animal phyla and markers with zero validation warnings:

| Organism & Phylum | Marker | Input Mode | Best-Fit BIC Model | Rooting |
| :--- | :--- | :--- | :--- | :--- |
| *Ancylostoma ceylanicum* (Nematoda) | Nuclear ITS | Raw `.ab1` (F + R) | IQ-TREE 2 MFP | Biological outgroup |
| *Allolobophora chlorotica* (Annelida) | Mitochondrial COI | Clean `.fasta` | `HKY+F+G4` | Midpoint (Run 1) |
| *Amolops* spp. (Amphibia) | Mitochondrial 16S rRNA | Clean `.fasta` + 8 refs | `TIM2+F+G4` | *Staurois parvus* (`AB259722`) |

---

## Quick-Start Protocol (Run 1 vs. Run 2 & Outgroup Guide)

1. **`.ab1` vs. `.fasta`:** Always upload **Forward `.ab1` (Box 1) + Reverse `.ab1` (Box 2)** for your own lab PCR samples from a sequencing company. Use `.fasta` only for already-trimmed consensus sequences or GenBank records.
2. **Chromatogram QC:** Keep Mott Phred Cutoff at `20` and inspect `01_chromatogram_qc.png` for crisp single peaks.
3. **Run 1 (Discovery):** Leave *Optional Published Reference Accessions* and *Optional Biological Outgroup* blank so BLASTn identifies your organism and exact gene length.
4. **Run 2 (Publication Figure):** Add 2–4 verified GenBank accessions from peer-reviewed papers + 1 **Biological Outgroup** following two golden rules:
   * **Same gene & similar length:** Match marker and length (`>=70–100%` of your query length).
   * **Close sister taxon:** Pick a closely related species in the same genus (or sister genus) so the tree is properly rooted without distorting ingroup branches.
5. **Publishing Your Results:** Upload `01_clean_query.fasta` to the [NCBI GenBank Submission Portal](https://submit.ncbi.nlm.nih.gov/about/genbank/) (requires `collection-date` and `geo_loc_name` metadata), and open `06_tree_for_MEGA_FigTree_iTOL.nwk` in [iTOL](https://itol.embl.de/) or [FigTree](https://github.com/rambaut/figtree/releases) for custom journal styling.

---

## Scope & Honest Limitations

* Designed specifically for **single-locus Sanger DNA barcoding** (e.g., ITS, COI, 16S, 18S, 28S)—not for NGS/WGS, metagenomics, or multi-gene supermatrix phylogenomics.
* Does not replace manual chromatogram inspection for heterozygous double peaks or noisy traces.
* Does not run Bayesian inference (MrBayes/BEAST) or automatically submit sequences to NCBI GenBank.

---

## Please Cite If You Use This Tool

This workbench is **100% free and open-source** for students, teachers, and researchers. If you use outputs from this app in a thesis, dissertation, poster, or peer-reviewed paper, **please cite it**—citations help keep free academic tools maintained and supported, and a quick GitHub **Star** or **Issue** with feedback/bug reports is always appreciated!

**Ready-to-paste Methods sentence:**
> Sanger sequence quality control, NCBI BLASTn reference curation, MAFFT multiple sequence alignment, ingroup population genetics, and IQ-TREE 2 Maximum-Likelihood phylogenetic inference were performed using Universal Sanger-to-Phylogeny Workbench v2.1 (Chettri, 2026; https://github.com/AbinayaChettri/sanger-phylogeny-workbench-NBBGC-TADONG-).

**Reference (APA / Journal style):**
> Chettri, A. (2026). *Universal Sanger-to-Phylogeny Workbench (v2.1)* [Computer software]. ORCID: https://orcid.org/0009-0000-9037-9699. GitHub: https://github.com/AbinayaChettri/sanger-phylogeny-workbench-NBBGC-TADONG- (Web application: https://nbbgc-tadong-phylogeny.streamlit.app).

**BibTeX:**
```bibtex
@software{chettri2026sanger_phylogeny_workbench,
  author       = {Chettri, Abinaya},
  title        = {Universal Sanger-to-Phylogeny Workbench},
  version      = {2.1},
  year         = {2026},
  url          = {https://github.com/AbinayaChettri/sanger-phylogeny-workbench-NBBGC-TADONG-},
  note         = {Live web application: https://nbbgc-tadong-phylogeny.streamlit.app; ORCID: 0009-0000-9037-9699}
}
```

You can also click **"Cite this repository"** in the right sidebar on GitHub (`CITATION.cff`).

---

## License

Released under the **MIT License** — Copyright (c) 2026 Abinaya Chettri.
