"""
### ======= Calculate read depth over the reference genome ======= ###
This snakemake file uses samtools to calculate the coverage for each position
in the wuhan reference genome (NC_045512.2) after various steps in the pipeline.

Author: Will Hannon 
"""

rule samtools_depth:
    """ 
    Calculate the depth over each position filtering by the phred base score. 
    """
    input: bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam")
    output: join(config['coverage_dir'], "{sample}", "{sample}.depth")
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


rule merge_depth:
    """ 
    Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: 
        expand(join(config['coverage_dir'], "{sample}", "{sample}.depth"), sample=['HL-2162'])
    output: 
        depth = join(config['coverage_dir'], "merged.depth.tsv"),
        header = temp(join(config['coverage_dir'], "merged.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Contig\tPOS\tDP\tRun"}}1' {output.header} > {output.depth}
        """


# Add a final rule to process the depth as a csv file
