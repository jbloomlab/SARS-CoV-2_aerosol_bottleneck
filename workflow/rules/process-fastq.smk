"""
### ======= Process the fastq files ======= ###
This snakemake file contains rules for processing raw fastqs.

# Author: Will Hannon 
"""

rule fetch_fastq:
    """ 
    Move the fastq files into the results directory and rename the files. 
    """
    output: R1 = join(config['fastq_dir'], "{sample}", "{sample}_R1.fastq.gz"),
            R2 = join(config['fastq_dir'], "{sample}", "{sample}_R2.fastq.gz")
    params: R1_path = lambda wildcards: sample_df.loc[(sample_df.Run == wildcards.sample)].R1.item(),
            R2_path = lambda wildcards: sample_df.loc[(sample_df.Run == wildcards.sample)].R2.item()
    shell: "cp {params.R1_path} {output.R1} && cp {params.R2_path} {output.R2}"


rule trim_adapters:
    """
    Use fastp to trim the adapters from the reads.

    1. Automatic adaptor trimming
    2. Low-qual base filtering (40% of bases w/ Phred >15) 
    3. Reporting by HTML and JSON.
    """
    input:  
        R1 = join(config['fastq_dir'], "{sample}", "{sample}_R1.fastq.gz"),
        R2 = join(config['fastq_dir'], "{sample}", "{sample}_R2.fastq.gz")
    output:
        R1 = join(config['trim_dir'], "{sample}", "{sample}_R1.fastq.gz"),
        R2 = join(config['trim_dir'], "{sample}", "{sample}_R2.fastq.gz"),
        html = join(config['qc_dir'], "{sample}", "{sample}.fastp.html"),
        json = join(config['qc_dir'], "{sample}", "{sample}.fastp.json")
    conda: '../envs/process-fastq.yml'
    shell:
        """ 
        fastp \
            -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
            --html {output.html} \
            --json {output.json} 
        """
