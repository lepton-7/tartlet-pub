# TaRTLEt Walkthrough

This tutorial for running the tool will use sample data that was part of the validation data set.

## Sample Data

All calls are relative to the repository root, function options and argument details may be accessed using the `--help` option.

Most of what is needed to run the following examples is shipped with the repository, except for the paired-end reads (due to size constraints). All the Sequence Run Archive (SRA) accessions used for this work are recorded in the `validation/rna_seq/acc_lists` directory with the `validation/rna_seq/download_SRAs.sh` convenience script used to download each dataset.

## Identify Riboswitches

TaRTLEt wraps around an infernal instance and comes packaged with rfam riboswitch covariance models. Calling:

```bash
tartlet-utils find-riboswitches -o validation/infernal/results --no-stats validation/genomes/*.fna
```

Runs infernal on the input genomes. Infernal output is a text table. To convert the results into a table that is used by TaRTLEt for further processing, call:

```bash
tartlet-utils make-table -i validation/infernal/results -o validation/tables/inf_results_test.csv
```

This table will be slightly different from that generated for the manuscript data (`validation/tables/inf_results.csv`) because of the `dataset` column. Values for this column were changed to be more concise than the tool outputs using a secondary script.

## Identify downstream ORFs

Information about open reading frames downstream of riboswitches is obtained by TaRTLEt by using prodigal. Calling:

```bash
tartlet-utils find-orfs -o validation/orfs/calls validation/genomes/*.fna
```

Runs prodigal on all the input genomes and:

```bash
tart-utils record-downstream --ledger validation/tables/inf_results_test.csv -i validation/orfs/calls
```

Adds ORFs downstream of riboswitches to the infernal results table. This information is used only to plot ORF locations in the locus plots.

## Extract Riboswitch information from input genome(s)

TaRTLEt's first step is to extract riboswitch sequences from input genomes using the infernal results table. The following steps of the pipeline will use the *B. subtilis* data as the example.

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

The `--dset` option instructs TaRTLEt to only consider the relevant microbe dataset from the cumulative results table. Next, the `bounds` command needs to be run to support proper location indexing during the next steps:

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

TaRTLEt uses HISAT2 for read alignment. The HISAT pipeline is mostly wrapped around by TaRTLEt. The alignment files are not included with the repository due to their size. Indexing:

```bash
tartlet-targeted index -i validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000 -p 44
```

The indexing input directory must be the directory with the `reference-gen` outputs. The alignment step is called directly through HISAT:

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

The SAM to sorted BAM conversion is wrapped by TaRTLEt:

```bash
tartlet-targeted convert-sam -i validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/alignment_final
```

## Process alignment results

```bash
tart-targeted parse-bam -i validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/alignment_final \
        -o validation/alignment/outputs/b_sub_168/plots/picks \
        --bounds-file validation/alignment/outputs/b_sub_168/switch_seqs_delta1000-1000/rowid_to_bounds.json \
        --allow-single-reads \
        --allow-soft-clips \
        --picks
```

## Generate TaRTLEt outputs

```bash
tart-targeted filter -i  validation/alignment/outputs/b_sub_168/plots/picks.tar.gz \
        -o validation/alignment/outputs/b_sub_168/plots/ \
        --ext-prop -0.3 1.0 \
        --conv \
        --min-cov-depth 15
```

The outputs of this are the `cluster_stats` and `peak_log` files placed in the output directory specified by the `-o` option.

## Plotting results (outside tool scope)

Currently, peak plots from the manuscript (Fig. 5) cannot be generated directly through the tool. The R script used to generate all the peak plots included in the manuscript and the supplement is located at `validation/plotting/peak_plotting.r`. The inputs to this script are the `cluster_stats` and `peak_log` files.
