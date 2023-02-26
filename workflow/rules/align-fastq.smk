"""
### ======= Align the reads to the wuhan reference genome ======= ###
This snakemake file uses BWA to align the trimmed reads to the wuhan reference genome.
These alignments are then sorted and indexed with samtools.

# Author: Will Hannon 
"""

rule get_reference:
    """ 
    Download the SARS-CoV-2 reference genome from the NCBI FTP site.
    """
    output: join(config['reference_dir'], 'SARS2.fa')
    params: ftp = config['ref']
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule get_annotations:
    """ 
    Download the SARS-CoV-2 gff annotations from the NCBI FTP site.
    """
    output: join(config['reference_dir'], 'SARS2.gff')
    params: ftp = config['gff']
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule bwa_index:
    """ 
    Index the genome with BWA before alignment with BWA. 
    """
    input: join(config['reference_dir'], 'SARS2.fa')
    output: join(config['reference_dir'], 'index', 'SARS2.fa')
    conda: '../envs/align-fastq.yml'
    params: algorithm="bwtsw"
    shell:
        """
        cp {input} {output}
        bwa index -a {params.algorithm} {output}
        """

        
rule bwa_align:
    """ 
    Perform short read alignment of the trimmed 
    viral readas with `bwa-mem`.
    
    Sort the aligned reads with samtools sort.
    """
    input: 
        reads = [join(config['trim_dir'], "{sample}", "{sample}_R1.fastq.gz"), 
                 join(config['trim_dir'], "{sample}", "{sample}_R2.fastq.gz")],
        genome = join(config['reference_dir'], 'index', 'SARS2.fa')
    output: bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
            bai = join(config['align_dir'], "{sample}", "{sample}.sorted.bam.bai")
    threads: config['threads']['max_cpu']
    conda: '../envs/align-fastq.yml'
    shell: 
        """
        bwa mem -t {threads} \
            {input.genome} \
            {input.reads} | \
            samtools view -bh | \
            samtools sort -o {output.bam} - 
        samtools index {output.bam}
        """

