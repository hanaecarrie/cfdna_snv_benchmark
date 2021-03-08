import re
import pysam
from tqdm import tqdm


def list_supporting_reads(bamfile_path, vcf_df):
    samfile = pysam.AlignmentFile(bamfile_path, "rb")

    # initiate list of reads to remove
    reads2remove = []
    log_dict = {"position": [], "type": [], "total_reads": [], 'supporting_reads': [],
                'normal_reads': [], 'alternative_reads': [], "problematic_reads": []}

    # iterate over positions
    for ci, mutation in tqdm(vcf_df.iterrows(), total=vcf_df.shape[0]):

        genotype = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        c = 0  # number of reads supporting the considered mutation
        t = 0  # total number of reads at that position
        p = 0  # number pf reads with issues
        n = 0  # number of reads supporting the reference genome
        a = 0  # number of reads with alternative nucleotide (not ref, not alt snp)
        mutation_type = None
        alt = []

        # iterate over reads that fall into the mutation position
        for read in samfile.fetch(str(mutation['#CHROM']), mutation['POS']-1, mutation['POS']):
            t += 1
            seq = read.query_alignment_sequence
            pos = (mutation['POS']-1) - read.reference_start
            cigar = read.cigarstring
            cond, cond1, cond2 = False, False, False

            # guess mutation type
            if ',' in mutation['ALT']:  # several possible variants
                if sum([len(muts) - len(mutation['REF']) == 0 for muts in mutation['ALT'].split(',')]):
                    mutation_type = 'SNV'
                elif sum([len(muts) - len(mutation['REF']) > 0 for muts in mutation['ALT'].split(',')]):
                    mutation_type = 'INS'
                elif sum([len(muts) - len(mutation['REF']) < 0 for muts in mutation['ALT'].split(',')]):
                    mutation_type = 'DEL'
                else:
                    raise ValueError('cannot identify mutation type: REF: ', mutation['REF'], 'ALT: ', mutation['ALT'])
            else:
                if len(mutation['ALT']) - len(mutation['REF']) == 0:
                    mutation_type = 'SNV'
                elif len(mutation['ALT']) - len(mutation['REF']) > 0:
                    mutation_type = 'INS'
                elif len(mutation['ALT']) - len(mutation['REF']) < 0:
                    mutation_type = 'DEL'
                else:
                    raise ValueError('cannot identify mutation type: REF: ', mutation['REF'], 'ALT: ', mutation['ALT'])

            ######## SNV ##########
            if mutation_type == 'SNV':
                # cond1 = True as cigar string is not informative of SNVs
                cond1 = True
                if cigar is not None:
                    cigar_pos = re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
                    cigar_states = re.split('[0-9]+', cigar)[1:]
                    cumul = 0
                    for i, cp in enumerate(cigar_pos):
                        if (cigar_states[i] != 'S') and (cumul <= pos):
                            cumul += -int(cp) if cigar_states[i] == 'D' else int(cp)
                            if cigar_states[i] == 'D':
                                pos += -int(cp)
                            elif cigar_states[i] == 'I':
                                pos += int(cp)
                genotype[seq[pos]] = genotype[seq[pos]]+1
                if ',' in mutation['ALT']:
                    for muts in mutation['ALT'].split(','):
                        if seq[pos] == muts:
                            cond2 = True
                else:
                    if seq[pos] == mutation['ALT']:
                        cond2 = True

            ######## INSERTION ##########
            elif mutation_type == 'INS':  # insertion
                # cond1 = cigar string indicates a deletion at this position
                if cigar is not None:
                    cigar_pos = re.split('M|I|D|N|S|H|P|=|X',cigar)[:-1]
                    cigar_states = re.split('[0-9]+', cigar)[1:]
                    if 'I' in cigar_states:
                        cigar_pos = [0 if (cigar_states[i] == 'S') else int(ci) for i, ci in enumerate(cigar_pos)]
                        cigar_pos = [-int(ci) if (cigar_states[i] == 'D') else int(ci) for i, ci in enumerate(cigar_pos)]
                        indexI = cigar_states.index('I')
                        if type(indexI) == list:
                            for idxI in indexI:
                                print(cigar_pos[:idxI], pos)
                                if sum(cigar_pos[:idxI]) == pos + 1:
                                    cond1 = True
                        else:
                            if sum(cigar_pos[:indexI]) == pos + 1:
                                cond1 = True
                if cond1:
                    # cond2 = nucleotide sequence comparison
                    if ',' in mutation['ALT']:
                        for muts in mutation['ALT'].split(','):
                            if seq[pos:pos + len(muts)] == muts:
                                cond2 = True
                    else:
                        cond2 = (seq[pos:pos+len(mutation['ALT'])] == mutation['ALT'])

            ######## DELETION ##########
            elif mutation_type == 'DEL':  # deletion
                # cond1 = cigar string indicates a deletion at this position
                if cigar is not None:
                    cigar_pos = re.split('M|I|D|N|S|H|P|=|X',cigar)[:-1]
                    cigar_states = re.split('[0-9]+',cigar)[1:]
                    if 'D' in cigar_states:
                        cigar_pos = [0 if (cigar_states[i] == 'S') else int(ci) for i, ci in enumerate(cigar_pos)]
                        cigar_pos = [-int(ci) if (cigar_states[i] == 'D') else int(ci) for i, ci in enumerate(cigar_pos)]
                        indexD = cigar_states.index('D')
                        if type(indexD) == list:
                            for idxD in indexD:
                                if sum(cigar_pos[:idxD]) == pos + 1:
                                    cond1 = True
                        else:
                            if sum(cigar_pos[:indexD]) == pos + 1:
                                cond1 = True
                if cond1:
                    # cond2 = nucleotide sequence comparison
                    if ',' in mutation['ALT']:
                        if seq[pos:pos + len(muts)] == muts:
                            cond2 = True
                    else:
                        cond2 = (seq[pos:pos+len(mutation['ALT'])] == mutation['ALT'])

            if cond1 and cond2:
                c += 1
                reads2remove.append(read.query_name)
            elif not cond1 and (seq[pos:pos+len(mutation['REF'])] == mutation['REF']):
                n += 1
            elif cigar is not None:
                a += 1
                alt.append(cond1)
                alt.append(seq[pos:pos+max(len(mutation['ALT']), len(mutation['REF']))])
            else:
                p += 1

        if mutation_type == 'SNV':
            print('SNV', mutation['POS'], 'REF:', mutation['REF'], 'ALT:', mutation['ALT'], genotype, '% VAF:', round(100*c/t, 2))
        else:
            print(mutation_type,  mutation['POS'], 'REF:', mutation['REF'], 'ALT:', mutation['ALT'], c, n, t, p, '% VAF:', round(100*c/t, 2), alt)
        log_dict["position"].append(mutation['POS'])
        log_dict["type"].append(mutation_type)
        log_dict["total_reads"].append(t)
        log_dict["supporting_reads"].append(c)
        log_dict["normal_reads"].append(n)
        log_dict["alternative_reads"].append(a)
        log_dict["problematic_reads"].append(p)

    samfile.close()

    return reads2remove, log_dict


if __name__ == '__main__':
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    bamfile_path = '../data/healthy_chr22_merged-ready.bam'
    vcf_df = pd.read_csv('../data/common_SNPs/dbsnp_df.csv')
    reads2remove, log_dict = list_supporting_reads(bamfile_path, vcf_df)
    print(len(reads2remove))
    log_pd = pd.DataFrame.from_dict(log_dict)
    quick_check = sum(log_pd['supporting_reads']) == len(reads2remove)
    print(quick_check)
    print('# reads to remove: ', len(reads2remove))
    print('% reads to remove: {:2f}%'.format(100*sum(log_pd['supporting_reads'])/sum(log_pd['total_reads'])))
    log_pd['vaf'] = log_pd['supporting_reads'] / log_pd['total_reads']
    log_pd['normal af'] = log_pd['normal_reads'] / log_pd['total_reads']
    log_pd['noisy af'] = log_pd['alternative_reads'] / log_pd['total_reads']
    plt.figure(figsize=(10, 5))
    plt.title('SNV')
    sns.histplot(data=log_pd[log_pd['type'] == 'SNV'][['vaf', 'normal af', 'noisy af']], bins=100,  stat="probability")
    plt.figure(figsize=(10, 5))
    plt.title('DEL')
    sns.histplot(data=log_pd[log_pd['type'] == 'DEL'][['vaf', 'normal af', 'noisy af']], bins=100,  stat="probability")
    plt.figure(figsize=(10, 5))
    plt.title('INS')
    sns.histplot(data=log_pd[log_pd['type'] == 'INS'][['vaf', 'normal af', 'noisy af']], bins=100,  stat="probability")
