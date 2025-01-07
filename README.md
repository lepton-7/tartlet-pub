# TaRTLEt Walkthrough

This is a step-by-step tutorial for running TaRTLEt on paired-end RNA-seq data to identify riboswitch loci showing evidence of condition-dependent transcription termination. In the tutorial, we'll use the *B. subtilis* subset of the validation dataset from our paper (Kshatriya and Bagby, "TaRTLEt: Transcriptionally-active Riboswitch Tracer Leveraging Edge deTection").

All paths given below are relative to the repository root. Function options and argument details may be accessed using the `--help` option.

## Download sample data

Due to size constraints, the paired-end reads and alignment files are not included in this repository. All the Sequence Run Archive (SRA) accessions used for this work are recorded in the `validation/rna_seq/acc_lists` directory. You can use the `validation/rna_seq/download_SRAs.sh` convenience script to download each dataset of interest. To download just the *B. subtilis* data, call

```bash
```

## Identify riboswitches

TaRTLEt wraps around an [Infernal](http://eddylab.org/infernal/) instance and comes packaged with [Rfam](https://rfam.org/) riboswitch covariance models. To identify riboswitch loci, call

```bash
tartlet-utils find-riboswitches -o validation/infernal/results --no-stats validation/genomes/*.fna
```

This runs Infernal on the input genomes, generating a text table. To convert the results into the format used by TaRTLEt for further processing, call

```bash
tartlet-utils make-table -i validation/infernal/results -o validation/tables/inf_results_test.csv
```

(Note that the format of the output's `dataset` column will differ slightly from that shown in `validation/tables/inf_results.csv`, where for concise presentation of the full validation dataset we used a secondary reformatting script.)

## Identify downstream ORFs

TaRTLEt's locus plots display the locations of open reading frames (ORFs) downstream of riboswitches. ORF annotations are also helpful downstream, for understanding the biological impacts of riboswitch activity. To identify these ORFs, TaRTLEt uses [Prodigal](https://github.com/hyattpd/Prodigal). To run Prodigal on all input genomes, call

```bash
tartlet-utils find-orfs -o validation/orfs/calls validation/genomes/*.fna
```

Then call

```bash
tart-utils record-downstream --ledger validation/tables/inf_results_test.csv -i validation/orfs/calls
```

to add downstream ORFs to the infernal results table.

## Extract riboswitch information from input genome(s)

TaRTLEt `reference-gen` uses the Infernal results table to extract riboswitch sequences from input genomes, retaining a window of user-definable size (`--pre-del` and `--post-del` parameters) around the riboswitch sequence. This extraction step improves computational efficiency in downstream read-mapping. We use the `--dset` option to specify a dataset to be extracted from the cumulative results table. To run riboswitch extraction, call

```bash
tartlet-targeted reference-gen \ 
        --ledger validation/tables/inf_results.csv \
        --out-dir validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000 \
        --genome validation/genomes \
        --dset b_sub_168 \
        --unify \
        --pre-del 1000 \
        --post-del 1000
```

Next, to support proper location indexing during the next steps, call `bounds`:

```bash
tartlet-targeted bounds \
        --ledger validation/tables/inf_results.csv \
        --out-dir validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000 \
        --genome validation/genomes \
        --dset b_sub_168 \
        --pre-del 1000 \
        --post-del 1000
```

## Run alignment

TaRTLEt uses HISAT2 for read alignment. First, index the riboswitch sequences by calling

```bash
tartlet-targeted index -i validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000 -p 44
```

(That is,  the `reference-gen` output directory is the input directory for indexing.) Then call HISAT to align reads:

```bash
hisat2 \
    -x validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/unified_index/unified_index \
    -1 <sequencing dir>/mate_1.fastq.gz \
    -2 <sequencing dir>/mate_2.fastq.gz \
    -S validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/alignment_final/unified/alignment.sam \
    -t \
    --no-unal \
    --score-min L,0,-0.4 \
    -p 40
```

Next, use TaRTLEt's `convert-sam` (a wrapper for HISAT2's SAM-to-sorted-BAM conversion function) to produce a sorted BAM file:

```bash
tartlet-targeted convert-sam -i validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/alignment_final
```

## Process alignment results

TaRTLEt `parse-bam` uses the sorted BAM files to produce its locus plots, compute the signed sum of read termini, and run convolution to identify candidate transcription-termination events. `parse-bam` also performs the first stage of significance testing, to identify which candidate events are above background biological and sequencing noise. This analysis is performed at the level of individual riboswitch loci, one experimental condition at a time. To process alignment results, call

```bash
tart-targeted parse-bam -i validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/alignment_final \
        -o validation/alignment/outputs/b_sub_168/plots/picks \
        --bounds-file validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/rowid_to_bounds.json \
        --allow-single-reads \
        --allow-soft-clips \
        --picks
```

This step generates many outputs:

- In the `fail/` subdirectory, one locus plot per riboswitch locus per condition for every case in which TaRTLEt detected no above-background transcription termination events. Locus plot filenames take the form `riboswitch_class#genome_accession_number#riboswitch_locus_start#riboswitch_locus_end#-#SRA_accession.sorted.bam.png`.
- In the `pass/` subdirectory, one locus plot per riboswitch locus per condition for every case in which TaRTLEt *did* detect above-background transcription termination event(s). Locus plot filenames take the form `riboswitch_class#genome_accession_number#riboswitch_locus_start#riboswitch_locus_end#-#SRA_accession.sorted.bam.png`.
- In the `picks/` subdirectory, extracted subsets of sorted BAM files. Loci are grouped by riboswitch family (e.g., `picks/TPP/`); each of these subdirectories contains a subdirectory for every SRA accession number (of the form `picks/riboswitch_family/SRA_accession.sorted.bam/`). The files in each subdirectory contain that condition's read-mapping results for each riboswitch locus in the family.

It can be helpful to inspect the `fail/` and `pass/` locus plots to convince yourself that these per-condition calls look reasonable for your dataset. The subsetted BAM files in `picks/` are not human-readable, but they're needed for the next step.

## Cluster peaks across conditions and test for condition-dependent transcription termination

TaRTLEt `filter` next compares read-mapping data across experimental conditions. The fractional coverage change across each peak in the convolution (per-locus, per-condition) is calculated, and peaks within the riboswitch region of interest (user-definable with `--ext-prop`) are clustered by x-position (per-locus, across-condition). Then for each cluster TaRTLEt asks (1) whether the mean change in coverage across the cluster is significantly more negative than expected by chance (Mann-Whitney one-tailed U-test) and (2) whether the variance in coverage change is significantly larger than expected by chance (Levene's test). You can set a minimum coverage depth threshold using `--min-cov-depth` and specify an output directory with `-o`. Here, we invoke `filter` as follows:

```bash
tart-targeted filter -i  validation/alignment/outputs/b_sub_168/plots/picks.tar.gz \
        -o validation/alignment/outputs/b_sub_168/plots/ \
        --ext-prop -0.3 1.0 \
        --conv \
        --min-cov-depth 15
```

This step generates two key outputs per riboswitch locus:

- `cluster_stats.csv` reports details on statistical tests of each cluster's fractional coverage change mean and variance.
- `peak_log.csv` reports peak-by-peak details, including the "noise sets" used for statistical comparisons.

## Plotting results (outside tool scope)

Currently, peak plots from the manuscript (Fig. 5) cannot be generated directly through the tool. The R script used to generate all the peak plots included in the manuscript and the supplement is located at `validation/plotting/peak_plotting.r`. The inputs to this script are the `cluster_stats` and `peak_log` files.
