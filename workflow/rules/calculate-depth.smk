"""
### ======= Calculate read depth over the reference genome ======= ###
This snakemake file uses samtools to calculate the coverage for each position
in the wuhan reference genome (NC_045512.2) after various steps in the pipeline.

Author: Will Hannon 
"""

rule samtools_untrimmed_depth:
    """ 
    Calculate the depth over each position filtering by the phred base score. 
    This is the depth before any primer trimming of the reads.
    """
    input: 
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam")
    output: 
        join(config['coverage_dir'], "{sample}", "{sample}.untrimmed.depth")
    params: 
        min_baseQ = config['variant_calling']['min_baseQ'],
        max_depth = config['variant_calling']['max_depth']
    conda: 
        '../envs/call-variants.yml'
    shell: 
        """
        samtools depth \
            -a \
            -q {params.min_baseQ} \
            -d {params.max_depth} \
            {input.bam} \
            > {output}

        sed -i "s/$/\t{wildcards.sample}/" {output} 
        """


rule merge_untrimmed_depth:
    """ 
    Merge the samtools depth tables for all of the accessions into a single file.
    This is the depth before any primer trimming of the reads.
    """
    input: 
        expand(join(config['coverage_dir'], "{sample}", "{sample}.untrimmed.depth"), sample=samples)
    output: 
        depth = join(config['coverage_dir'], "merged.untrimmed.depth.tsv"),
        header = temp(join(config['coverage_dir'], "merged.untrimmed.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Contig\tPOS\tDP\tRun"}}1' {output.header} > {output.depth}
        """


rule samtools_trimmed_depth:
    """ 
    Calculate the depth over each position filtering by the phred base score. 

    This is after trimming the reads to remove the primer sequences.
    """
    input: 
        bam = join(config['align_dir'], "{sample}", "{sample}.masked.sorted.bam")
    output: 
        join(config['coverage_dir'], "{sample}", "{sample}.trimmed.depth")
    params: 
        min_baseQ = config['variant_calling']['min_baseQ'],
        max_depth = config['variant_calling']['max_depth']
    conda: 
        '../envs/call-variants.yml'
    shell: 
        """
        samtools depth \
            -a \
            -q {params.min_baseQ} \
            -d {params.max_depth} \
            {input.bam} \
            > {output}

        sed -i "s/$/\t{wildcards.sample}/" {output} 
        """


rule merge_trimmed_depth:
    """ 
    Merge the samtools depth tables for all of the accessions into a single file.

    This is after trimming the reads to remove the primer sequences.
    """
    input: 
        expand(join(config['coverage_dir'], "{sample}", "{sample}.trimmed.depth"), sample=samples)
    output: 
        depth = join(config['coverage_dir'], "merged.trimmed.depth.tsv"),
        header = temp(join(config['coverage_dir'], "merged.trimmed.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Contig\tPOS\tDP\tRun"}}1' {output.header} > {output.depth}
        """

