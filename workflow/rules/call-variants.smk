"""
### ======= Call intrahost viral variants ======= ###
This snakemake file uses multiple approaches to call intrahost
viral variants from the aligned reads. It attemps to standardize
hueristic filters when possible .

Author: Will Hannon 
"""

rule get_varscan:
    """ 
    Download Varscan into tools/ directory. Fetch from github. 
    """
    output: join(config['tool_dir'], "VarScan.v2.4.0.jar")
    params: http=config['varscan']
    shell: "wget -O {output} {params.http}"


rule varscan_calling:
    """ 
    SNP calling with Varscan 2 using an mpileup file generated with samtools.
    """
    input: 
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa'),
        varscan = join(config['tool_dir'], "VarScan.v2.4.0.jar")
    output: 
        variants = join(config['variant_dir'], "{sample}", "{sample}.varscan.vcf")
    params: 
        max_depth = config['variant_calling']['max_depth'],
        min_baseQ = config['variant_calling']['min_baseQ'],
        min_mapQ = config['variant_calling']['min_mapQ'],
        min_coverage = config['variant_calling']['min_coverage'],
        min_frequency = config['variant_calling']['min_frequency'],
        min_obsv = config['variant_calling']['min_obsv'],
        strand_filter = config['variant_calling']['strand_filter']
    conda: 
        '../envs/call-variants.yml'
    shell: 
        """
        samtools mpileup \
            -aa \
            -A \
            -d {params.max_depth} \
            -B \
            -Q 0 \
            -q {params.min_mapQ} \
            -f {input.genome} \
            {input.bam} | \
        java -jar {input.varscan} mpileup2snp - \
            --output-vcf 1 \
            --min-coverage {params.min_coverage} \
            --min-reads2 {params.min_obsv} \
            --min-avg-qual {params.min_baseQ} \
            --strand-filter {params.strand_filter} \
            --min-var-freq {params.min_frequency} > {output.variants}
        """


rule ivar_calling:
    """ 
    SNP calling with iVar using an mpileup file generated with samtools.
    """
    input: 
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa'),
        gff = join(config['reference_dir'], 'SARS2.gff')
    output: 
        variants = join(config['variant_dir'], "{sample}", "{sample}.ivar.tsv")
    params: 
        max_depth = config['variant_calling']['max_depth'],
        min_baseQ = config['variant_calling']['min_baseQ'],
        min_mapQ = config['variant_calling']['min_mapQ'],
        min_coverage = config['variant_calling']['min_coverage'],
        min_frequency = config['variant_calling']['min_frequency']
    conda: 
        '../envs/call-variants.yml'
    shell: 
        """
        samtools mpileup \
            -aa \
            -A \
            -d {params.max_depth} \
            -B \
            -Q 0 \
            -q {params.min_mapQ} \
            -f {input.genome} \
            {input.bam} | \
        ivar variants \
            -p {output.variants} \
            -q {params.min_baseQ} \
            -t {params.min_frequency} \
            -m {params.min_coverage} \
            -r {input.genome} \
            -g {input.gff}
        """


rule lofeq_calling: 
    """
    Call vaiants with lofreq using a probabilistic model. 
    """
    input: 
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa')
    output: 
        variants = join(config['variant_dir'], "{sample}", "{sample}.lofreq.vcf")
    params: 
        max_depth = config['variant_calling']['max_depth'],
        min_baseQ = config['variant_calling']['min_baseQ'],
        min_mapQ = config['variant_calling']['min_mapQ'],
        min_coverage = config['variant_calling']['min_coverage']
    threads: config['threads']['max_cpu']
    conda: 
        '../envs/call-variants.yml'
    shell: 
        """       
        lofreq call-parallel --pp-threads {threads} \
            -f {input.genome} \
            -d {params.max_depth} \
            -q {params.min_baseQ} \
            -Q {params.min_baseQ} \
            -m {params.min_mapQ} \
            -C {params.min_coverage} \
            {input.bam} \
            -o {output}
        """


rule get_SnpEff:
    """ 
    Download and build SnpEff with annotations from gtf file. 
    """
    output: 
        join(config['tool_dir'], 'snpEff/snpEff.jar')
    params:
        installdir = config['tool_dir'],
        installpath = config['snpeff'],
        tmp = join(config['tool_dir'], 'snpEff_latest_core.zip'),
        datadir = join(config['tool_dir'], 'snpEff/data')
    shell:
        """
        # Download SnpEff
        wget {params.installpath} -O {params.tmp}

        # Unzip and install
        unzip {params.tmp} \
            -d {params.installdir} \
            && rm -rf {params.tmp}

        # Set path variable to snpEff dir
        snpeff={output}

        # Make the data directory to install viral genomes
        mkdir -p {params.datadir}
        """


rule build_SnpEff:
    """ 
    Build the annotation repository for the genome used to call variants. 
    """
    input: 
        snpeff = join(config['tool_dir'], 'snpEff/snpEff.jar'),
        ref = join(config['reference_dir'], 'SARS2.fa'),
        gff = join(config['reference_dir'], 'SARS2.gff')
    output: 
        virusdir = directory(join(config['tool_dir'], 'snpEff/data/SARS2'))
    params:
        config=join(config['tool_dir'], 'snpEff/snpEff.config'),
        genome="SARS2"
    conda: 
        '../envs/call-variants.yml'
    shell:
        """
        # Make the virus specific dir
        mkdir -p {output.virusdir}

        # Move the genome fasta file
        cp {input.ref} {output.virusdir}/sequences.fa

        # Move the gtf file
        cp {input.gff} {output.virusdir}/genes.gff

        # Add the genome build to the end of the file
        echo "{params.genome}.genome: {params.genome}" >> {params.config}

        # Build the new database
        java -Xss100M -Xmx8g -jar {input.snpeff} build -gff3 -v {params.genome}
        """


rule annotate_vcf:
    """
    Use the program SnpEff to annotate the effect of 
    mutations using the custom virus genome as reference.

    This tool is set up in the download_tools.smk file.
    """
    input:
        virusdir = join(config['tool_dir'], 'snpEff/data/SARS2'),
        vcf = join(config['variant_dir'], "{sample}", "{sample}.{caller}.vcf")
    output: 
        join(config['variant_dir'], "{sample}", "{sample}.{caller}.ann.vcf")
    conda: 
        '../envs/call-variants.yml'
    params:
        snpEff = join(config['tool_dir'], 'snpEff/snpEff.jar'),
        config = join(config['tool_dir'], "snpEff/snpEff.config"),
        genome = "SARS2"
    shell:
        """
        java -jar {params.snpEff} -c {params.config} -noStats \
            -v {params.genome} {input.vcf} > \
            {output}
        """


rule vcf_to_table:
    """
    Convert the VCF files to tables for easy data
    analysis in R or Python.
    """
    input: 
        join(config['variant_dir'], "{sample}", "{sample}.{caller}.ann.vcf")
    output: 
        join(config['variant_dir'], "{sample}", "{sample}.{caller}.tsv")
    conda: 
        '../envs/call-variants.yml'
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F QUAL -F REF -F ALT \
            -F DP -F AF -F FILTER -GF DP \
            -GF RD -GF FREQ -GF SDP -GF AD -F ANN \
            -O {output}
        """
