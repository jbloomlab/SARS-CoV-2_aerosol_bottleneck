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
    SNP calling with Varscan v2 using an mpileup file generated with samtools.
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
        strand_filter = config['variant_calling']['strand_filter'],
        min_pval = config['variant_calling']['min_pval']
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
            --p-value {params.min_pval} \
            --strand-filter {params.strand_filter} \
            --min-var-freq {params.min_frequency} > {output.variants}
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
        min_coverage = config['variant_calling']['min_coverage'],
        min_pval = config['variant_calling']['min_pval']
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
            -a {params.min_pval} \
            {input.bam} \
            -o {output}
        """


rule pysam_calling: 
    """
    Custom variant calling script using pysam htslib interface. 
    """
    input: 
        bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
        genome = join(config['reference_dir'], 'index', 'SARS2.fa'),
        gff = join(config['reference_dir'], 'SARS2.gff')
    output:
        variants = join(config['variant_dir'], "{sample}", "{sample}.pysam.formatted.tsv")
    params: 
        max_depth = config['variant_calling']['max_depth'],
        min_baseQ = config['variant_calling']['min_baseQ'],
        min_coverage = config['variant_calling']['min_coverage'],
        min_frequency = config['variant_calling']['min_frequency']
    conda:
        '../envs/call-variants.yml'
    shell:
        """
        python workflow/scripts/pysam_variant_caller.py \
            -b {input.bam} \
            -r {input.genome} \
            -g {input.gff} \
            -d {params.max_depth} \
            -q {params.min_baseQ} \
            -c {params.min_coverage} \
            -f {params.min_frequency} \
            -o {output.variants}
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


rule format_ivar:
    """
    Format the ivar output to be compatible with the other variant callers.
    """
    input: 
        ivar = join(config['variant_dir'], "{sample}", "{sample}.ivar.tsv")
    output: 
        formatted = join(config['variant_dir'], "{sample}", "{sample}.ivar.formatted.tsv")
    params:
        min_pval = config['variant_calling']['min_pval']
    run:
        # Read in the iVar variants 
        df = pd.read_csv(str(input.ivar), sep="\t")
        # Remove the failing variants
        df = df[df['PVAL'] >= params.min_pval]
        # Rename the columns to keep
        df.rename(columns={
                    'REF_RV': 'RDR',
                    'ALT_RV': 'ADR', 
                    'REF_QUAL': 'RBQ',
                    'ALT_QUAL': 'ABQ',
                    'TOTAL_DP': 'DP',
                    'ALT_FREQ': 'AF'
                    }, inplace=True)
        # Add in missing columns
        df['RDF'] = df['REF_DP'] - df['RDR']
        df['ADF'] = df['ALT_DP'] - df['ADR']
        df['GENE'] = df['GFF_FEATURE'].str.split(':').str[0]
        df['CALLER'] = 'ivar'
        # Filter out rows with indels
        df = df[df['ALT'].str.len() <= 1]
        # Drop unecessary columns 
        df = df.drop(columns=[
            'REGION', 'REF_DP', 'ALT_DP', 'PVAL',
            'PASS', 'REF_CODON', 'ALT_CODON', 'GFF_FEATURE'])
        # Drop duplicate rows
        df.drop_duplicates(inplace=True)
        # Write the output
        df.to_csv(str(output.formatted), sep="\t", index=False)


rule get_SnpEff:
    """ 
    Download and build SnpEff with annotations from gff file. 
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


rule vcf_to_tsv:
    """
    Convert the VCF files to tables for downstream analysis.
    """
    input: 
        join(config['variant_dir'], "{sample}", "{sample}.{caller}.ann.vcf")
    output: 
        join(config['variant_dir'], "{sample}", "{sample}.{caller}.formatted.tsv")
    conda: 
        '../envs/call-variants.yml'
    shell:
        """
        python workflow/scripts/vcf_to_tsv.py \
            -i {input} \
            -c {wildcards.caller} \
            -a \
            -o {output}
        """


rule merge_variants:
    """
    Aggregate the variant calls from each caller into a single table.
    Remove variants that fail universal filters.
    Remove variants in masked sites (excluded sites table).
    """
    input: 
        variants = [
            join(config['variant_dir'], "{sample}", "{sample}.lofreq.formatted.tsv"),
            join(config['variant_dir'], "{sample}", "{sample}.varscan.formatted.tsv"),
            join(config['variant_dir'], "{sample}", "{sample}.ivar.formatted.tsv"),
            join(config['variant_dir'], "{sample}", "{sample}.pysam.formatted.tsv")
        ]
    output:
        join(config['variant_dir'], "{sample}", "{sample}.variants.tsv")
    params:
        exclude = config['excluded'],
        min_coverage = config['variant_calling']['min_coverage'],
        min_frequency = config['variant_calling']['min_frequency'],
        min_obsv = config['variant_calling']['min_obsv']
    run:
        # Read in each variant caller from input list
        dfs = [pd.read_csv(f, sep="\t") for f in input.variants]
        # Concatenate the dataframes
        df = pd.concat(dfs, ignore_index=True)
        # Add the sample name as a column
        df['Run'] = wildcards.sample
        # Apply final filters that are universal to all callers
        df = df[df['AF'] >= params.min_frequency]
        df = df[df['DP'] >= params.min_coverage]
        df[df['ADF'] + df['ADR'] >= params.min_obsv]
        # Read in the excluded sites table (from the Lauring lab)
        exclude = pd.read_csv(params.exclude, sep="\t")
        # Remove variants where the position is in the excluded sites table
        df = df[~df['POS'].isin(exclude['POS'])]
        # Write the output
        df.to_csv(str(output), sep="\t", index=False)

# Add a rule to aggregate all variants into a single csv file