# Reproducible sheep 10x BCR V(D)J integration

Code-only repository for a two-stage integration of sheep 10x 5′ BCR sequencing. **No custom Cell Ranger reference, IgMAT reference, primers, FASTQs, or subject-level sequence data are distributed.** Supply your own validated reference and primer sequences. All Python scripts use the Python 3.9+ standard library.

## Study design

| Sample | Sequencing 1 | Sequencing 2 IgG | Sequencing 2 IgM |
|---|---|---|---|
| sample1 | PUR1 | 1G | 1M |
| sample2 | PUR2 | 2G | 2M |

Barcodes are comparable **within each sample** because libraries derive from the same initial 10x barcoded preparation. A barcode in sample1 is **not** the same cell as the identical barcode in sample2. Run labels do not determine isotype: constant-region `c_gene` annotations do.

## Quick start

```bash
# 1. Provide FASTQ locations and your own Cell Ranger-compatible sheep VDJ reference.
bash scripts/run_cellranger_vdj.sh --reference /your/validated/sheep_vdj_reference --manifest config/cellranger_libraries.csv --outdir results/cellranger
# The above previews commands. Add --execute only when paths are verified.

# 2. Place the six filtered_contig_annotations.csv files under data/ as named in config/libraries.csv.
python3 scripts/integrate_vdj.py --manifest config/libraries.csv --outdir results/integration

# 3. Primer mapping, using your actual primers in CSV columns primer_id,sequence:
python3 scripts/primer_map_fasta.py --input sequences.fasta --primers config/my_primers.csv --output results/fasta_primer_hits.tsv --max-mismatches 1
python3 scripts/primer_map_fastq.py --input reads.fastq.gz --primers config/my_primers.csv --output results/fastq_primer_hits.tsv --max-mismatches 1
```

`config/primers.example.csv` contains an **illustrative, fictional** primer and must not be used as a study primer. Primer matching supports IUPAC ambiguity codes, both orientations, zero-based coordinates and optional substitutions; it does not perform alignment with gaps, adapter trimming, or paired-end read assembly. Report the primer list, matching orientation, allowed mismatches, and upstream preprocessing in the manuscript.

## Outputs

For each of `seq1`, `seq2`, `merged`:
- `*_cells.tsv`: one row per unique `(sample, barcode)` and exclusive chain category.
- `*_unique_chains.tsv`: distinct productive chains with run provenance, V/J genes, all seven FR/CDR AA and NT regions, reconstructed variable-region AA and NT, translation QC.
- `*_strict_1H1L_pairs.tsv`: **one row per unambiguous cell with exactly one distinct IGH and one distinct IGK/IGL**, heavy/light fields side by side.
- `*_all_candidate_pairs.tsv`: includes **candidate** H/L combinations in multichain cells; not all represent verified biological pairs.

Also:
- `summary_observed_vs_unique.tsv`: samplewise library-level category observations (sum before integration) versus post-merge unique barcode categories, including IgA/IgD/IgE/Other.
- `seq2_barcode_overlap.tsv`: G/M shared and run-exclusive barcodes.
- `barcode_category_changes.tsv`: classification changes from sequencing 1 and 2 to final merge.

## Rules and caveats

1. Input rows must have `is_cell=true`, `productive=true`, and chain `IGH`, `IGK`, or `IGL`. Cell Ranger filtered annotations are not equivalent to all raw reads.
2. Deduplication within `(sample, barcode)` uses `(chain, reconstructed FR1–FR4 variable-region nucleotide sequence, c_gene)` **only when translation QC passes**. For incomplete/failed reconstructions, the conservative fallback includes the contig ID, preventing potentially unsafe collapse. Therefore unresolved contigs may inflate multichain counts and require manual QC.
3. Cell categories are mutually exclusive. Exactly one heavy and one light is a strict 1H1L pair; multiple distinct chains are classified as multichain. Isotype is determined by the prefix of `c_gene`, with unknown/missing annotations classified as `Other`.
4. The reconstructed nucleotide sequence is **FR1 through FR4**, not the full original contig, leader, or constant region. A single extra terminal nucleotide in Cell Ranger `fwr4_nt` is removed when necessary to align with annotated `fwr4` AA. Translation is checked against annotated amino acids. `translation_qc` flags mismatches and missing regions.
5. Library-level counts and merged unique category counts are **not category-wise subsets**: new heavy/light pairings and multichain discoveries can change classifications. Overall observations are not independent biological cells.
6. Different constant-gene annotations on an otherwise identical heavy variable sequence are retained as separate chain records. This is conservative for potential isotype ambiguity and may require biological review.
7. Cell Ranger commands assume the installed Cell Ranger version supports `cellranger vdj` with `--chain=IG`; record and cite the actual version used. Do not claim an IgMAT processing stage unless its exact code/reference is released or described separately.
8. Before submission, run against the original six CSVs, review translation-QC flags and compare counts to the archived analysis. Record software versions and checksums in the manuscript repository release.

## Reproducibility and data availability

The scripts are reusable for any comparable 10x V(D)J dataset after editing the manifest. They do **not** redistribute the study's unpublished custom reference or the original sequencing data. Publish the source FASTQ accession and a versioned GitHub release/Zenodo DOI separately, subject to repository and journal policy. A GitHub repository URL alone is mutable; a tagged release and archived DOI are preferable.

## Suggested repository citation

`Authors. (YEAR). Sheep 10x BCR repertoire processing and barcode-level integration (v1.0.0). GitHub repository and archived DOI [to be added].`

## License

MIT license; see LICENSE.
