"""
The purpose of this script is to call variants from a BAM file using pysam.

Author: Will Hannon
"""
from collections import defaultdict
import argparse
import pysam
import pandas as pd
from Bio import SeqIO


def check_read(read):
    """
    Helper function to decide what reads should
    be keep when parsing alignment file.

    Parameters
    ----------
    read : pysam.AlignedSegment
        read from an aligned BAM file

    Returns
    -------
    bool
        True/False if read should be included

    """
    # Exclude Quality Failures
    if read.is_qcfail:
        return False
    # Exclude Secondary Mappings
    elif read.is_secondary:
        return False
    # Exclude Unmapped Reads
    elif read.is_unmapped:
        return False
    # Exclude Marked Duplicates
    elif read.is_duplicate:
        return False
    else:
        return True


def call_variants(bampath, reference, contig, max_Depth=100000, min_BQ=25, min_Freq=0.005, min_Cov=1, callback=check_read):
    """
    Call variants from a pileup file using pysam.


    Parameters
    ----------
    bampath : str
        Path to an aligned, sorted, and indexed BAM file
    reference : str
        The reference sequence from a parsed FASTA file
    contig : str
        The chromosome name from a parsed FASTA file
    min_BQ : int
        The minimum base quality for a read to be considered at a position
    min_Cov : int
        The minimum depth for a position to be considered
    callback : function
        A function that return True if a read should be considered


    Returns
    -------
    pd.DataFrame
        A table for variants identified in a given sample


    """
    # Store the rows
    rows = []

    # Read in the indexed BAM file
    with pysam.AlignmentFile(bampath, "rb") as bamfile:
        # Iterate over each column (a.k.a. position) in the pileup
        for pileupcolumn in bamfile.pileup(contig, stepper='nofilter', min_base_quality=min_BQ, max_depth=max_Depth):

            # Site-specific information
            POS = pileupcolumn.reference_pos + 1
            REF = reference[pileupcolumn.reference_pos]
            DP = 0
            RDF = 0
            RDR = 0
            RBQ = []
            ALT = {base: {'ADF': 0,
                          'ADR': 0,
                          'ABQ': []
                          } for base in 'ATCG' if base != REF}

            # Log a progress message
            if POS % 1000 == 0:
                print(
                    f"Called variants for {POS} of {len(reference)} total bases...")

            # Iterate over each read in the pileup
            for pileupread in pileupcolumn.pileups:
                # Skip reads that fail the QC requirements or contain indels
                if pileupread.is_refskip or pileupread.is_del or not callback(pileupread.alignment):
                    continue

                # Increment Total Depth
                DP += 1

                # Collect Read-specific information
                readpos = pileupread.query_position
                base = pileupread.alignment.query_sequence[readpos]
                qual = pileupread.alignment.query_qualities[readpos]
                orientation = 'reverse' if pileupread.alignment.is_reverse else 'forward'

                # Check if the site is the reference base
                if base == REF:
                    RBQ.append(qual)
                    if orientation == "forward":
                        RDF += 1
                    else:
                        RDR += 1
                # Or, fill in SNPs at this position
                elif base != REF:
                    ALT[base]['ABQ'].append(qual)
                    if orientation == "forward":
                        ALT[base]['ADF'] += 1
                    else:
                        ALT[base]['ADR'] += 1

            # Collect all of the information for this positions
            for SNP, SNP_stats in ALT.items():
                # If there is a SNP at this position, add it to the table
                if SNP_stats['ADF'] + SNP_stats['ADR'] > 0:
                    # Calculate the mean alternative base quality
                    MEAN_ABQ = round(
                        sum(SNP_stats['ABQ'])/len(SNP_stats['ABQ']), 2)
                    # Calculate the mean reference base quality (it can be zero)
                    if RDF + RDR == 0:
                        MEAN_RBQ = 0
                    else:
                        MEAN_RBQ = round(sum(RBQ)/len(RBQ), 2)
                    AF = (SNP_stats['ADF'] + SNP_stats['ADR']) / DP
                    row = [POS, REF, SNP, round(
                        AF, 4), DP, RDF, RDR, SNP_stats['ADF'], SNP_stats['ADR'], MEAN_RBQ, MEAN_ABQ]
                    rows.append(row)

    # Convert the rows into a dataframe
    variant_df = pd.DataFrame(rows, columns=[
                              "POS", "REF", "ALT", "AF", "DP", "RDF", "RDR", "ADF", "ADR", "RBQ", "ABQ"])
    variant_df['CALLER'] = 'pysam'

    # Apply any remaining filters
    variant_df = variant_df.loc[variant_df["AF"] >= min_Freq]
    variant_df = variant_df.loc[variant_df["DP"] >= min_Cov]

    return variant_df.reset_index(drop=True)


def translate(codon):
    """
    Translate a three letter DNA string into 
    a one letter amino acid code. 

    Parameters
    ----------
    codon : str
        three letter DNA sequence

    Returns
    -------
    str
        one letter amino acid code

    Raises
    ------
    AssertionError
        error if codon sequence is invalid

    """

    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    assert codon in table.keys(), "Not a valid codon sequence."

    return table[codon]


def parse_gff_file(filename):
    """
    Generic parsing function for GFF3 format files.

    Parameters
    ----------
    filename : str
        Path to GFF3 format annotation file

    Returns
    -------
    dict
        Dictionary of parsed features and header

    """

    data = {
        'header': {},
        'features': []
    }

    with open(filename) as f:
        print(f"Parsing {filename}...\n")

        for line in f:
            line = line.strip()

            if line.startswith('###'):
                print("End of file.")
                break
            elif line.startswith('##'):
                key, value = line[2:].split(' ', 1)
                data['header'][key] = value
            elif line.startswith('#'):
                continue
            else:
                fields = line.split('\t')
                feature = {
                    'contig': fields[0],
                    'source': fields[1],
                    'type': fields[2],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'score': fields[5],
                    'strand': fields[6],
                    'phase': fields[7],
                    'attributes': {}
                }

                # Parse the attributes field
                for attribute in fields[8].split(';'):
                    if not attribute:
                        continue
                    key, value = attribute.split('=')
                    feature['attributes'][key] = value

                data['features'].append(feature)

    return data


def annotate_coding_change(position, alt, gff, reference):
    """
    Given a mutation and alternative base, annotate the coding 
    effect if there is one.

    Parameters
    ----------
    position : int
        1-indexed integer position of mutation in reference genome

    alt : str
        Alternative base at this position

    gff : dict
        Parsed dictionary of GFF annotations 

    reference : str
        Reference sequence from fasta file


    Returns
    -------
    list
        List of [REF_AA, PROT_POS, ALT_AA, GENE]

    """

    # Extract the coding regions from the GFF file
    CDS = defaultdict(list)
    for feature in gff['features']:
        if feature['type'] != 'CDS':
            continue
        # Append the gene name, start and stop (1-indexed)
        gene = feature['attributes']['gene']
        start = feature['start']
        stop = feature['end']
        CDS[gene].append([start, stop])

    # Annotate any possible coding changes
    coding_changes = []
    for gene, regions in CDS.items():
        # Note that one mutation can have multiple coding effects in a gene (overlapping RFs and frameshifts)
        for [start, stop] in regions:
            if position >= start and position <= stop:
                # Position of the SNP in the gene (0-index)
                gene_pos = position - start
                # Position of the SNP in the codon (0-index)
                codon_pos = gene_pos % 3
                # Position of the ALT_AA in the protein (0-index)
                prot_pos = (gene_pos // 3)
                # Gene sequence
                gene_seq = str(reference)[start-1:stop]
                # List of the codons in the gene
                codons = [gene_seq[i:i+3] for i in range(0, len(gene_seq), 3)]
                # Reference Codon
                ref_codon = codons[prot_pos]
                # Alternative Codon
                alt_codon = "".join(
                    alt if (i == codon_pos) else base for i, base in enumerate(ref_codon))
                # Reference AA
                ref_aa = translate(ref_codon)
                # Alternative AA
                alt_aa = translate(alt_codon)

                # Append the coding changes
                coding_changes.append([ref_aa, prot_pos + 1, alt_aa, gene])

    # Join multiple annotations if they exist
    if coding_changes:
        return [":".join(str(x) for x in y) if len(set(y)) > 1 else y[0] for y in zip(*coding_changes)]
    else:
        return [None, None, None, None]


def main():

    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        description='Call variants from a BAM file using Pysam.')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='Path to input BAM file sorted and indexed')
    parser.add_argument('-r', '--reference', type=str, required=True,
                        help='Path to a reference genome in FASTA format')
    parser.add_argument('-g', '--gff', type=str, required=True,
                        help='Path to annotation file in GFF3 format')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to output TSV file')
    parser.add_argument('-q', '--min-baseQ', type=int, default=25,
                        help='Minimum base quality to consider a read')
    parser.add_argument('-d', '--max-depth', type=int, default=100000,
                        help='Maximum depth of reads per position')
    parser.add_argument('-f', '--min-freq', type=float, default=0.005,
                        help='Minimum frequency to call a report a variant')
    parser.add_argument('-c', '--min-coverage', type=int, default=1,
                        help='Minimum coverage to include a position')
    # Parse the arguments
    args = parser.parse_args()

    # Parse the GFF file
    gff = parse_gff_file(args.gff)

    # Parse the reference genome
    refrecord = SeqIO.read(args.reference, "fasta")
    contig = refrecord.name
    reference = refrecord.seq

    # Call variants
    variant_df = call_variants(args.bam,
                               reference,
                               contig,
                               min_BQ=args.min_baseQ,
                               min_Freq=args.min_freq,
                               min_Cov=args.min_coverage,
                               max_Depth=args.max_depth,
                               callback=check_read)

    # Annotate the coding changes
    variant_df[['REF_AA', 'POS_AA', 'ALT_AA', 'GENE']] = variant_df.apply(
        lambda row: pd.Series(annotate_coding_change(row['POS'], row['ALT'], gff, reference)), axis=1)

    # Write the output to a TSV file
    variant_df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
