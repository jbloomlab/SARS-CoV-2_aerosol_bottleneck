"""
The purpose of this script is to convert a VCF file to a TSV file.

Author: Will Hannon
"""
import re
import argparse
import vcf
import pandas as pd


def one_letter_format(aa):
    """
    Convert 3 letter amino acid format to 1 letter amino acid format

    Parameters
    ----------
    aa: str
        Amino acid in three letter format

    Returns
    -------
    str
        Amino acid in one letter format

    """
    amino_acids = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
        'Asx': 'B', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
        'Glx': 'Z', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
        'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
        'Tyr': 'Y', 'Val': 'V', '*': '*', 'Ter': '*'
    }
    return amino_acids[aa]


def vcf_to_table(vcf_path, caller, annotations=False):
    """
    Read a vcf file created with either lofreq or varscan
    into a pandas dataframe.

    Parameters
    ----------
    vcf_path : str
        Path to a VCFv4.0/VCFv4.1 format file from lofreq or varscan
    caller : ['lofreq', 'varscan']
        The variant calling program was used to generate the vcf file
    annotations : bool
        Shoud annotations made with SnpEFF be parsed from the file

    Returns
    -------
    pd.DataFrame
        A pandas dataframe containing the key fields

    """

    # Make sure that caller is a viable option
    if caller not in {'lofreq', 'varscan'}:
        raise ValueError(
            f"{caller} is not one of the supported options: ['lofreq', 'varscan']")

    # Make an empty list to store vcf information
    vcf_info = []

    # Read in the vcf file from provided path
    with open(vcf_path, 'r') as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)

        # Iterate over each record in the file
        for record in vcf_reader:

            # Get the fields common to all formats
            POS = record.POS
            REF = record.REF
            ALT = record.ALT[0]

            # Get infomration from lofreq formatted file
            if caller == "lofreq":
                DP = record.INFO['DP']
                AF = record.INFO['AF']
                RDF, RDR, ADF, ADR = record.INFO['DP4']
                RBQ = None
                ABQ = None

            # Get the information from a varscan formatted file
            if caller == "varscan":
                # There will be only a single sample per record
                for sample in record.samples:
                    DP = sample['DP']
                    AF = float(sample['FREQ'].strip('%')) / 100
                    RDF = sample['RDF']
                    RDR = sample['RDR']
                    ADF = sample['ADF']
                    ADR = sample['ADR']
                    RBQ = sample['RBQ']
                    ABQ = sample['ABQ']

            # Get the annotations if they're present
            if annotations:
                ANN = record.INFO['ANN'][0].split("|")
                GENE = ANN[3]
                match = re.match(
                    r'^p\.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})$', ANN[10])
                if match:
                    REF_AA = one_letter_format(match.group(1))
                    POS_AA = match.group(2)
                    ALT_AA = one_letter_format(match.group(3))
                else:
                    REF_AA = None
                    POS_AA = None
                    ALT_AA = None
                # Add the fields from this record
                vcf_info.append([POS, REF, ALT, AF, DP, RDF, RDR, ADF,
                                ADR, RBQ, ABQ, GENE, REF_AA, POS_AA, ALT_AA])
            # If there aren't annotations, return the full
            else:
                # Add the fields from this record
                vcf_info.append(
                    [POS, REF, ALT, AF, DP, RDF, RDR, ADF, ADR, RBQ, ABQ])

    # Write the info out to a dataframe
    if annotations:
        vcf_df = pd.DataFrame(vcf_info,
                              columns=["POS", "REF", "ALT", "AF", "DP", "RDF", "RDR", "ADF", "ADR", "RBQ", "ABQ", "GENE", "REF_AA", "POS_AA", "ALT_AA"])
        vcf_df['CALLER'] = caller
    else:
        vcf_df = pd.DataFrame(vcf_info,
                              columns=["POS", "REF", "ALT", "AF", "DP", "RDF", "RDR", "ADF", "ADR", "RBQ", "ABQ"])
        vcf_df['CALLER'] = caller

    return vcf_df


def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        description='Convert a lofreq or varscan fromat VCF file into a TSV file.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to input VCF file in lofreq or varscan format')
    parser.add_argument('-c', '--caller', type=str, required=True,
                        help='Which variant caller was used to generate the VCF file: (lofreq or varscan))')
    parser.add_argument('-a', '--add_annotations', action='store_true',
                        help='Add annotations to output file')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to output TSV file')
    # Parse the arguments
    args = parser.parse_args()

    # Convert the vcf file to a tsv file
    vcf_df = vcf_to_table(args.input, args.caller, args.add_annotations)

    # Write the dataframe to a tsv file
    vcf_df.to_csv(args.output, sep='\t', index=False)


if __name__ == "__main__":
    main()
