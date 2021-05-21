"""
Filtering mutated positions present in VCF file from a BAM file
"""
import io
import sys
import pysam
import pandas as pd


INPUT_VCF = str(sys.argv[1])
INPUT_BAM = str(sys.argv[2])
OUTPUT_BAM = str(sys.argv[3])

print(INPUT_VCF)
print(INPUT_BAM)
print(OUTPUT_BAM)


def get_all_positions(vcf_file):
    """
    Generates a sorted list of positions from input SMURF file

    Arguments:
        vcf_file {str} -- path to VCF output file

    Return:
        {list(tup)} -- tup[0] = chromosome
                       tup[1] = start_pos
                       tup[2] = stop_pos
    """
    with open(vcf_file, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    vcf_posis = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    vcf_posis = vcf_posis[['CHROM', "POS"]].set_index('CHROM')
    vcf_posis['POS'] = vcf_posis['POS'].astype(str)
    vcf_posis = vcf_posis.groupby('CHROM').agg({'POS': ', '.join }).to_dict()['POS']
    for k, v in vcf_posis.items():
        vcf_posis[k] = list(map(int, v.split(", ")))
    
    return vcf_posis

def generate_pileup_posis(vcf_posis, thresh):
    """
    generate regions to look for when piling up in cancer file

    Arguments:
         vcf_posis {{[list]}} -- dict of lists containing chromosome as key and positions sorted
    """
    pileups = []
    for key in vcf_posis:
        for i, _ in enumerate(vcf_posis[key]):
            if i != 0:
                if (vcf_posis[key][i] - 150) - (vcf_posis[key][i-1] + 150) > thresh:
                    pileups.append((key, vcf_posis[key][i-1] + 150 + int(thresh/2), vcf_posis[key][i] - 150 - int(thresh/2)))

    return pileups

def run_pileups(samfile, output, pileups):
    with pysam.AlignmentFile(output, "wb", header=samfile.header) as outf:
        for p in pileups:
            reads = samfile.fetch(p[0], p[1], p[2])

            for read in reads:
                outf.write(read)

def main():
    """
    Main workflow
    """
    samfile = pysam.AlignmentFile(INPUT_BAM, "rb")
    positions = get_all_positions(INPUT_VCF)
    pileups = generate_pileup_posis(positions, 10)
    run_pileups(samfile, OUTPUT_BAM, pileups)
    samfile.close()


if __name__ == "__main__":
    main()
