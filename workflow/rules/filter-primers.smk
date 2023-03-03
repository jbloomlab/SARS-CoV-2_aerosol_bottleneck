"""
### ======= Remove atrifacts due to amplicon sequencing ======= ###
This snakemake file uses iVar to remove amplicon sequencing artifacts from
the aligned reads. It uses the primer bed file to identify the amplicon regions.

This workflow is based on the following resources:
    - https://doi.org/10.1186/s13059-018-1618-7
    - https://doi.org/10.1038/s41467-023-36001-5
    - https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

Author: Will Hannon 
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
        primer_pairs = pd.read_csv(input.primers)
        with open(output.fasta, "w") as f:
            for i, row in primer_pairs.iterrows():
                f.write(f">{row['Name']}" + "\n")
                f.write(f"{row['Sequence']}" + "\n")


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


rule test_primer_trimming:
    """
    Debugging rule to test primer trimming.
    """
    input:
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam")
    output:
        bam = join(config['align_dir'], "{sample}", "{sample}.subsampled.sorted.bam"),
        bai = join(config['align_dir'], "{sample}", "{sample}.subsampled.sorted.bam.bai")
    conda:
        "../envs/filter-primers.yml"
    shell:
        """
        samtools view -s 0.01 -b {input.bam} | \
            samtools sort -o {output.bam} && \
            samtools index {output.bam}
        """


rule trim_primers:
    """
    Use iVar to trim the primer sequences from the aligned reads.
    """
    input:
        # bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
        bed = join(config['reference_dir'], "primers", "primers.bed"),
        pairs = config['primer_pairs'],
        bam = join(config['align_dir'], "{sample}", "{sample}.subsampled.sorted.bam")
    output:
        bam = join(config['align_dir'], "{sample}", "{sample}.trimmed.bam")
    log:
        join(config['qc_dir'], "{sample}", "{sample}.ivar.trim.log")
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        ivar trim \
            -i {input.bam} \
            -b {input.bed} \
            -f {input.pairs} \
            -e \
            -p {output.bam} \
            &> {log}
        """


rule sort_trimmed_bams:
    """
    Sort and index the primer trimmed bam files with samtools.
    """
    input:
        bam = join(config['align_dir'], "{sample}", "{sample}.trimmed.bam")
    output:
        bam = join(config['align_dir'], "{sample}", "{sample}.trimmed.sorted.bam"),
        bai = join(config['align_dir'], "{sample}", "{sample}.trimmed.sorted.bam.bai")
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
        bam = join(config['align_dir'], "{sample}", "{sample}.trimmed.sorted.bam"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa')
    output:
        consensus = join(config['consensus_dir'], "{sample}", "{sample}.consensus.fa")
    params:
        max_depth = config['consensus_calling']['max_depth'],
        min_depth = config['consensus_calling']['min_depth'],
        min_baseQ = config['consensus_calling']['min_baseQ'],
        freq_threshold = config['consensus_calling']['freq_threshold']
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
            -q {params.min_baseQ} \
            -t {params.freq_threshold} \
            -m {params.min_depth}
        """


rule index_consensus:
    """
    Index the sample specific consensus sequence with bwa and samtools.
    """
    input:
        join(config['consensus_dir'], "{sample}", "{sample}.consensus.fa")
    output:
        join(config['consensus_dir'], "{sample}", "index", "{sample}.consensus.fa")
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
        genome = join(config['consensus_dir'], "{sample}", "index", "{sample}.consensus.fa")
    output:
        bam = join(config['consensus_dir'], "{sample}", "primers", "{sample}.primers.bam"),
        bed = join(config['consensus_dir'], "{sample}", "primers", "{sample}.primers.bed")
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
        bam = join(config['consensus_dir'], "{sample}", "primers", "{sample}.primers.bam"),
        genome = join(config['consensus_dir'], "{sample}", "index", "{sample}.consensus.fa")
    output:
        mismatches = join(config['consensus_dir'], "{sample}", "primers", "{sample}.mismatches.tsv")
    params:
        max_depth = config['primer_calling']['max_depth'],
        min_baseQ = config['primer_calling']['min_baseQ'],
        min_mapQ = config['primer_calling']['min_mapQ'],
        min_frequency = config['primer_calling']['min_frequency']
    conda: 
        "../envs/filter-primers.yml"
    shell:
        """
        samtools mpileup \
            -A \
            -d {params.max_depth} \
            --reference {input.genome} \
            -Q {params.min_baseQ} \
            -q {params.min_mapQ} \
            {input.bam} | \
        ivar variants \
            -p {output.mismatches} \
            -t {params.min_frequency}
        """


rule mask_primer_mismatches:
    """
    Use iVar to mask primer mismatches.
    """
    input:
        mismatches = join(config['consensus_dir'], "{sample}", "primers", "{sample}.mismatches.tsv"),
        bed = join(config['consensus_dir'], "{sample}", "primers", "{sample}.primers.bed"),
        pairs = config['primer_pairs']
    output:
        mask = join(config['consensus_dir'], "{sample}", "primers", "{sample}.masked.primers.txt")
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
        bam = join(config['align_dir'], "{sample}", "{sample}.trimmed.sorted.bam"),
        bed = join(config['consensus_dir'], "{sample}", "primers", "{sample}.primers.bed"),
        mask = join(config['consensus_dir'], "{sample}", "primers", "{sample}.masked.primers.txt")
    output:
        bam = join(config['align_dir'], "{sample}", "{sample}.masked.bam"),
        # bam = join(config['align_dir'], "{sample}", "{sample}.masked.sorted.bam"),
        # bai = join(config['align_dir'], "{sample}", "{sample}.masked.sorted.bai")
    conda: 
        "../envs/filter-primers.yml"
    log:
       join(config['qc_dir'], "{sample}", "{sample}.ivar.removereads.log")
    shell:
        """
        ivar removereads \
            -i {input.bam} \
            -p {output.bam} \
            -t {input.mask} \
            -b {input.bed} \
            &> {log}
        """


rule sort_masked_bams:
    """
    Sort and index the primer mismatch masked bam files with samtools.
    """
    input:
        bam = join(config['align_dir'], "{sample}", "{sample}.masked.bam")
    output:
        bam = join(config['align_dir'], "{sample}", "{sample}.masked.sorted.bam"),
        bai = join(config['align_dir'], "{sample}", "{sample}.masked.sorted.bam.bai")
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