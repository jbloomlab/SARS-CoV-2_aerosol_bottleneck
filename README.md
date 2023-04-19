# Determining the infulence of aerosol size on SARS-CoV-2 transmission dynamics

This analysis was conducted by William Hannon, Jesse Bloom, and Hui-Ling Yen. The manuscript associated with this analysis is avaliable [here]().

## Methods

### Sequencing Read Processing

We processed the illumina sequencing reads using a Snakemake pipeline that is available on GitHub [here](https://github.com/jbloomlab/SARS-CoV-2_aerosol_bottleneck). Briefly, the program [fastp](https://github.com/OpenGene/fastp) was used to trim reads of adaptor sequences and filter reads with an abundance of low quality bases (> 40% of bases with a phred score < 15). Following read trimming and filtering, the FASTQs were aligned to the Wuhan-1/2019 reference (NC_045512.2) using [BWA mem](https://bio-bwa.sourceforge.net/bwa.shtml). To avoid PCR artifacts that could bias variant calling, we used [iVar](https://andersen-lab.github.io/ivar/html/) to trim primer sequences from the aligned sequencing reads. Finally, we assessed the quality of aligned and trimmed BAM files by checking the read coverage with the program [Samtools](http://www.htslib.org/).

### Variant Calling and Filtering

To identify intrahost single nucleotide variants (SNVs), we used a [custom variant calling tool written in python](https://github.com/jbloomlab/SARS-CoV-2_aerosol_bottleneck/blob/main/workflow/scripts/pysam_variant_caller.py). Briefly, we used the python/samtools interface [pysam](https://github.com/pysam-developers/pysam) to iterate through each position in the aligned genome and keep track of the read depth, number of reference bases, number of alternative bases, read orientation, and base quality. The output of the tool is a row for each SNV along with its frequency, strand bias, number of observations, and coding effect. In addition to our own custom approach, we used three other variant calling programs – [lofreq](https://csb5.github.io/lofreq/), [iVar](https://andersen-lab.github.io/ivar/html/), and [varscan2](https://varscan.sourceforge.net/) – to identify intrahost SNVs. We annotated the coding effect of SNVs identified with these tools using [SnpEff](http://pcingola.github.io/SnpEff/).

Whenever possible, we standardized the filters between all four variant calling approaches. If a given filter could not be applied within a program, we applied the filter post hoc in R. When calling SNVs, we capped the maximum depth over each position at 100,000 reads and included only bases with a phred base quality > 25. We removed SNVs that were identified in fewer than 2% of reads or SNVs that had fewer than 200 reads covering a position. Overall, the concordance between approaches was high with > 85% of SNVs identified by all four programs (_Supplementary Figure X_).

Each sample apart from the isolate used to inoculate the animals was sequenced in technical replicates. The SNV allele frequencies had good concordance between replicates (_Supplementary Figure Y_). To create a single set of SNVs per sample, we kept only SNVs that were identified in both technical replicates and took the allele frequency from the replicate with the highest read depth. Following this, we finalized a single set of variants between variant calling approaches by keeping the SNVs identified by our custom variant calling program that were also identified by at least one other variant calling program.

### Bottleneck Estimation

To avoid overestimating the bottleneck size due to false-positive SNVs that arise due to sequencing errors, PCR amplification, or computational artifacts, we applied a more stringent set of filters. We removed SNVs that were observed on fewer than 25 reads, SNVs that were identified on more than 90% of reads of the same strand orientation, SNVs that were present in primer sequences, and SNVs that were identified in only a single batch of sequencing runs.

We then defined donor-contact transmission pairs based on the design of the experiment. In the case where a single cage had two inoculated donor hamsters connected to a cage with two naive contact hamsters, we considered all four pairwise combinations of donors and contacts. In the case where a single inoculated hamster was used as the donor for a single naive hamster, this was the only transmission pair.

To calculate the transmission bottleneck between each individual donor-contact pair, we adapted [an implementation](https://github.com/weissmanlab/BB_bottleneck/blob/master/Bottleneck_size_estimation_approx.r) of the approximate beta-binomial model that was initially defined by [Sobel Leonard, et. al, 2017](https://journals.asm.org/doi/10.1128/JVI.00171-17).

### Code Availability

All of the code used to run the analysis are available on GitHub at https://github.com/jbloomlab/SARS-CoV-2_aerosol_bottleneck. The code is also archived on Zenodo at [DOI:XXXXXXXXXXXXX]().

### Data Availability

The raw sequencing reads are available on the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) under the BioProject accession number [PRJNAXXXXXX]().
