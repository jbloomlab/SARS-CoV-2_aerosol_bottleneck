"""
### ======= Remove atrifacts due to amplicon sequencing ======= ###
This snakemake file uses iVar to remove amplicon sequencing artifacts from
the aligned reads. It uses the primer bed file to identify the amplicon regions.

This workflow is based on the following resources:
    - https://doi.org/10.1186/s13059-018-1618-7
    - https://doi.org/10.1038/s41467-023-36001-5
    - https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

# Author: Will Hannon 
"""

rule make_primer_fasta:
    """
    Create a fasta file of the primer sequences.
    """
    input:
        primers = config['primers']
    output:
        fasta = join(config['reference_dir'], "primers", "primers.fasta")
    run:
        primer_positions = pd.read_csv(input.primers)
        with open(output.fasta, "w") as f:
            for i, row in primer_positions.iterrows():
                f.write(f">{row['start']}_{row['end']}_{row['chromosome'][0:-2]}_{row['direction']}" + "\n")
                f.write(f"{row['sequence']}" + "\n")


rule make_primer_pairs:
    """
    TODO: Get the primer names for Hui-Ling's HKU Primer Set.
    """
    input:
        primers = config['primers']
    output:
        tsv = join(config['reference_dir'], "primers", "pairs.tsv")
    run:
        primer_positions = pd.read_csv(input.primers)


rule make_primer_bed:
    """
    Create a bed file of the primer sequences.
    """
    input:
        fasta = join(config['reference_dir'], "primers", "primers.fasta"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa')
    output:
        bam = temp(join(config['reference_dir'], "primers", "primers.bam")),
        bed = join(config['reference_dir'], "primers", "primers.bed")
    threads:
        config['threads']['max_cpu']
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        bwa mem -t {threads} \
            -k 5 \
            -T 16 \
            {input.genome} \
            {input.fasta} | \
            samtools view -b -F 4 > \
            {output.bam}
        
        bedtools bamtobed -i {output.bam} > {output.bed}
        """


rule trim_primers:
    """
    Use iVar to trim the primer sequences from the aligned reads.
    """
    input:
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
        bed = join(config['reference_dir'], "primers", "primers.bed")
    output:
        bam = join(config['primertrim_dir'], "{sample}", "{sample}.trimmed.bam")
    log:
        join(config['qc_dir'], "{sample}", "{sample}.ivar.log")
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        ivar trim \
            -i {input.bam} \
            -b {input.bed} \
            -p {output.bam} \
            &> {log}
        """


rule sort_trimmed_bams:
    """
    Sort and index the primer trimmed bam files with samtools.
    """
    input:
        bam = join(config['primertrim_dir'], "{sample}", "{sample}.trimmed.bam")
    output:
        bam = join(config['primertrim_dir'], "{sample}", "{sample}.trimmed.sorted.bam"),
        bai = join(config['primertrim_dir'], "{sample}", "{sample}.trimmed.sorted.bam.bai")
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        samtools sort \
            -o {output.bam} \
            {input.bam} && \
        samtools index \
            {output.bam}
        """


rule get_sample_consensus:
    """
    Use samtools and iVar to get a consensus sequence from the 
    aligned and primer-trimmed reads.
    """
    input:
        bam = join(config['primertrim_dir'], "{sample}", "{sample}.trimmed.sorted.bam"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa')
    output:
        consensus = join(config['ivar_dir'], "{sample}", "consensus", "{sample}.consensus.fa")
    params:
        max_depth = 100000,
        min_depth = 10,
        min_qual = 0,
        freq_threshold = 0
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
    	samtools mpileup \
            -a \
            -A \
            -d {params.max_depth} \
            -Q 0 \
            --reference {input.genome} \
            {input.bam} | \
        ivar consensus \
            -p {output.consensus} \
            -n N \
            -q {params.min_qual} \
            -t {params.freq_threshold} \
            -m {params.min_depth}
        """


rule index_consensus:
    """
    Index the sample specific consensus sequence with bwa and samtools.
    """
    input:
        join(config['ivar_dir'], "{sample}", "consensus", "{sample}.consensus.fa")
    output:
        join(config['ivar_dir'], "{sample}", "consensus", "index", "{sample}.consensus.fa")
    params: 
        algorithm="bwtsw"
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        cp {input} {output}
        bwa index -a {params.algorithm} {output}
        samtools faidx {output}
        """


rule make_consensus_primer_bed:
    """
    Create a bed file of the primer sequences for each sample consensus sequence.
    """
    input:
        fasta = join(config['reference_dir'], "primers", "primers.fasta"),
        genome = join(config['ivar_dir'], "{sample}", "consensus", "index", "{sample}.consensus.fa")
    output:
        bam = join(config['ivar_dir'], "{sample}", "primers", "{sample}.primers.bam"),
        bed = join(config['ivar_dir'], "{sample}", "primers", "{sample}.primers.bed")
    threads:
        config['threads']['max_cpu']
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        bwa mem -t {threads} \
            -k 5 \
            -T 16 \
            {input.genome} \
            {input.fasta} | \
            samtools view -b -F 4 | \
            samtools sort -o {output.bam}
        
        bedtools bamtobed -i {output.bam} > {output.bed}
        """


rule call_primer_mismatches:
    """
    Use iVar to identify sample specific primer mismatches.
    """
    input:
        bam = join(config['ivar_dir'], "{sample}", "primers", "{sample}.primers.bam"),
        consensus = join(config['ivar_dir'], "{sample}", "consensus", "index", "{sample}.consensus.fa")
    output:
        mismatches = join(config['ivar_dir'], "{sample}", "mismatches", "{sample}.mismatches.tsv")
    params:
        max_depth = 100000,
        min_base_qual = 30,
        min_map_qual = 20,
        freq_threshold = 0.02
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        samtools mpileup \
            -aa \
            -A \
            -d {params.max_depth} \
            --reference {input.consensus} \
            -Q {params.min_base_qual} \
            -q {params.min_map_qual} \
            -F 0 \
            {input.bam} | \
        ivar variants \
            -p {output.mismatches} \
            -t {params.freq_threshold}
        """


rule mask_primer_mismatches:
    """
    Use iVar to mask primer mismatches.
    """
    input:
        mismatches = join(config['ivar_dir'], "{sample}", "mismatches", "{sample}.mismatches.tsv"),
        bed = join(config['ivar_dir'], "{sample}", "primers", "{sample}.primers.bed"),
        pairs = join(config['reference_dir'], "primers", "pairs.tsv")
    output:
        mask = join(config['ivar_dir'], "{sample}", "mask", "{sample}.masked.primers.txt")
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        ivar getmasked \
            -i {input.mismatches} \
            -b {input.bed} \
            -f {input.pairs} \
            -p {output.mask} 
        """


rule remove_primer_mismatches:
    """
    Use iVar to remove reads with primer mismatches from trimmed bam file.
    """
    input:
        bam = join(config['primertrim_dir'], "{sample}", "{sample}.trimmed.sorted.bam"),
        bed = join(config['ivar_dir'], "{sample}", "primers", "{sample}.primers.bed"),
        mask = join(config['ivar_dir'], "{sample}", "mask", "{sample}.masked.primers.txt")
    output:
        bam = join(config['ivar_dir'], "{sample}", "{sample}.mismatched.trimmed.sorted.bam"),
        bai = join(config['ivar_dir'], "{sample}", "{sample}.mismatched.trimmed.sorted.bai")
    params:
        remove_sites_1 = "data/ivar_output/removed/{sample}_1.masked",
        temp_1 = "data/ivar_output/removed/{sample}_1.tmp"
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        ivar removereads \
            -i {input.bam} \
            -p {params.remove_sites_1} \
            -t {input.mask} \
            -b {input.bed}  

        samtools sort \
            -T {params.temp_1} \
            -o {output.bam} \
            {params.remove_sites_1}.bam && \
        samtools index {output.bam}
        """
