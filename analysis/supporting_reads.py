import re
import pysam
import pandas as pd
from tqdm.notebook import tqdm


def list_reads_to_remove(bamfile_path, common_snps_df, patient_snps_df, max_vaf=0.1, verbose=0):
    samfile = pysam.AlignmentFile(bamfile_path, "rb")

    # initiate list of reads to remove
    reads2remove = []
    log_dict = {'#CHROM': [], "POS": [], "type": [], 'ID': [], 'REF': [], 'ALT': [],
                "total_reads": [], 'supporting_reads': [], 'normal_reads': [], 'alternative_reads': [], "problematic_reads": [],
                "is patient snp": [], "patient snp alt": [], "patient snp vaf": []}

    # list positions in patient
    listpatientsnps = get_listpatientsnps(patient_snps_df)

    # iterate over known snps positions
    for ci, mutation in tqdm(common_snps_df.iterrows(), total=common_snps_df.shape[0]):

        genotype = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}  # genotype if SNV
        c = 0  # number of reads supporting the considered mutation
        t = 0  # total number of reads at that position
        p = 0  # number pf reads with issues
        n = 0  # number of reads supporting the reference genome
        a = 0  # number of reads with alternative nucleotide (not ref, not alt snp)
        mutation_type = get_mutation_type(mutation['REF'], mutation['ALT'])  # determine mutation type
        reads2remove_tmp = []  # temporary list of reads to remove at that locus, check VAF before adding to final list

        # iterate over reads that fall into the mutation position
        # /!\ the vcf file is 1-index based but pysam is 0-index based
        for read in samfile.fetch(str(mutation['#CHROM']), mutation['POS']-1, mutation['POS']):
            t += 1
            seq = read.query_alignment_sequence
            pos = (mutation['POS']-1) - read.reference_start
            cigar = read.cigarstring
            cond = False

            ######## SNV ##########
            if mutation_type == 'SNV':
                # cond = sequence matching
                # the mutation position needs to account for the reads indels
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
                            cond = True
                else:
                    if seq[pos] == mutation['ALT']:
                        cond = True

            ######## INSERTION / DELETION ##########
            elif (mutation_type == 'INS') or (mutation_type == 'DEL'):
                # cond = cigar string indicates a deletion at this position
                if cigar is not None:
                    cigar_pos = re.split('M|I|D|N|S|H|P|=|X',cigar)[:-1]
                    cigar_states = re.split('[0-9]+', cigar)[1:]
                    if mutation_type[0] in cigar_states:
                        cigar_pos = [0 if (cigar_states[i] == 'S') else int(ci) for i, ci in enumerate(cigar_pos)]
                        cigar_pos = [-int(ci) if (cigar_states[i] == 'D') else int(ci) for i, ci in enumerate(cigar_pos)]
                        indexINDEL = cigar_states.index(mutation_type[0])
                        if type(indexINDEL) == list:
                            for idxINDEL in indexINDEL:
                                print(cigar_pos[:idxINDEL], pos)
                                if sum(cigar_pos[:idxINDEL]) == pos + 1:
                                    cond = True
                        else:
                            if sum(cigar_pos[:indexINDEL]) == pos + 1:
                                cond = True
            if cond:  # read supporting known variant
                c += 1
                reads2remove_tmp.append(read.query_name)
            elif not cond and (seq[pos:pos+len(mutation['REF'])] == mutation['REF']):  # normal read
                n += 1
            elif cigar is not None:  # read supporting an unknown variant
                a += 1
            else:  # problematic read that does not have a cigar string, unmapped
                p += 1
        if verbose > -1:
            if mutation_type == 'SNV':
                print('SNV', mutation['POS'], 'REF:', mutation['REF'], 'ALT:', mutation['ALT'], genotype, c, n, t, p, '% VAF:', round(100*c/t, 2))
            else:
                print(mutation_type,  mutation['POS'], 'REF:', mutation['REF'], 'ALT:', mutation['ALT'], c, n, t, p, '% VAF:', round(100*c/t, 2))
        if t > 0:
            if round(c/t, 2) <= max_vaf:  # if VAF <= 10% (default), few reads supporting variants
                if mutation['ID'] not in listpatientsnps:  # snps found in healthies is not in patient's snps
                    for r in reads2remove_tmp:
                        reads2remove.append(r)  # remove rare mutated reads
                else:
                    # if same variant => do not remove, if other variant => remove
                    if mutation['ALT'] != patient_snps_df[patient_snps_df['ID'].str.contains(mutation['ID'])]['ALT'].values[0]:
                        for r in reads2remove_tmp:
                            reads2remove.append(r)  # remove rare alternative variant

        log_dict["#CHROM"].append(mutation['#CHROM'])
        log_dict["POS"].append(mutation['POS'])
        log_dict["type"].append(mutation_type)
        log_dict["total_reads"].append(t)
        log_dict["supporting_reads"].append(c)
        log_dict["normal_reads"].append(n)
        log_dict["alternative_reads"].append(a)
        log_dict["problematic_reads"].append(p)
        log_dict["ID"].append(mutation['ID'])
        log_dict["REF"].append(mutation['REF'])
        log_dict["ALT"].append(mutation['ALT'])
        log_dict["is patient snp"].append(mutation['ID'] in list(patient_snps_df['ID']))
        log_dict["patient snp alt"].append(patient_snps_df[patient_snps_df['ID'] == mutation['ID']]['ALT'].values[0]
                                           if mutation['ID'] in list(patient_snps_df['ID']) else '')
        log_dict["patient snp vaf"].append(patient_snps_df[patient_snps_df['ID'] == mutation['ID']]['VAF'].values[0]
                                           if mutation['ID'] in list(patient_snps_df['ID']) else '')

    log_pd = pd.DataFrame.from_dict(log_dict)
    print('# reads to remove: ', len(reads2remove))
    print('% reads to remove: {:2f}%'.format(100*sum(log_pd['supporting_reads'])/sum(log_pd['total_reads'])))
    log_pd['vaf'] = log_pd['supporting_reads'] / log_pd['total_reads']
    log_pd['normal af'] = log_pd['normal_reads'] / log_pd['total_reads']
    log_pd['noisy af'] = log_pd['alternative_reads'] / log_pd['total_reads']

    samfile.close()

    return reads2remove, log_pd


def prepare_bamsurgeon_inputs(patient_snps_df, log_pd, max_vaf=0.1):

    listpatientsnps = get_listpatientsnps(patient_snps_df)
    bamsurgeon_snv_dict = {'chr': [], 'pos_start': [], 'pos_end': [], 'vaf': [], 'alt': []}  # 1-index based
    bamsurgeon_indel_dict = {'chr': [], 'pos_start': [], 'pos_end': [], 'vaf': [], 'type': [], 'alt': []}  # 0-index based

    # iterate over mutated loci in the healthies with high VAF (above 10%)
    for ci, mutation in tqdm(log_pd[log_pd['vaf'] > max_vaf].iterrows(), total=log_pd[log_pd['vaf'] > max_vaf].shape[0]):
        patient_snp = patient_snps_df[patient_snps_df['ID'].str.contains(mutation['ID'])].squeeze()
        condA = mutation['ID'] not in listpatientsnps
        condB = mutation['ID'] in listpatientsnps and mutation['ALT'] in patient_snp['ALT']
        if condA or condB:
            # if not in patient's SNPs or a different variant from the patient's snps
            # bamsurgeon REF genome with VAF = 1
            bamsurgeon_snv_dict, bamsurgeon_indel_dict = add_mutation_bamsurgeon_dict(
                bamsurgeon_snv_dict, bamsurgeon_indel_dict, mutation, vaf=1, alt=mutation['REF'])

    # iterate through patients SNPs
    for pi, patient_snp in tqdm(patient_snps_df.iterrows(), total=patient_snps_df.shape[0]):
        # bamsurgeon patient genotype with VAF of the patient
        bamsurgeon_snv_dict, bamsurgeon_indel_dict = add_mutation_bamsurgeon_dict(
            bamsurgeon_snv_dict, bamsurgeon_indel_dict, patient_snp, patient_snp['VAF'], patient_snp['ALT'])

    bamsurgeon_snv_pd = pd.DataFrame.from_dict(bamsurgeon_snv_dict)
    bamsurgeon_indel_pd = pd.DataFrame.from_dict(bamsurgeon_indel_dict)

    return bamsurgeon_snv_pd, bamsurgeon_indel_pd


def add_mutation_bamsurgeon_dict(bamsurgeon_snv_dict, bamsurgeon_indel_dict, mutation, vaf, alt):
    mutation_type = get_mutation_type(mutation['REF'], alt)
    if ',' not in alt:
        if mutation_type == 'SNV':
            bamsurgeon_snv_dict['chr'].append(str(mutation['#CHROM']))
            bamsurgeon_snv_dict['pos_start'].append(mutation['POS'])
            bamsurgeon_snv_dict['pos_end'].append(mutation['POS'])
            bamsurgeon_snv_dict['vaf'].append(vaf)
            bamsurgeon_snv_dict['alt'].append(alt[0])
        else:  # mutation is insertion or deletion
            if mutation_type == 'DEL':
                len_mut = len(mutation['REF']) - len(alt)  # important for BamSurgeon
            else:  # mutation_type == 'INS'
                len_mut = len(alt) - len(mutation['REF'])  # not important for BamSurgeon
            bamsurgeon_indel_dict['chr'].append(str(mutation['#CHROM']))
            bamsurgeon_indel_dict['pos_start'].append(mutation['POS']-1)
            bamsurgeon_indel_dict['pos_end'].append(mutation['POS']+len_mut-1)
            bamsurgeon_indel_dict['vaf'].append(vaf)
            bamsurgeon_indel_dict['type'].append(mutation_type)
            bamsurgeon_indel_dict['alt'].append(alt if mutation_type == 'INS' else '')
    else:  # multiple variants
        for mi, mut in enumerate(alt.split(',')):
            if mutation_type[mi] == 'SNV':
                bamsurgeon_snv_dict['chr'].append(str(mutation['#CHROM']))
                bamsurgeon_snv_dict['pos_start'].append(mutation['POS'])
                bamsurgeon_snv_dict['pos_end'].append(mutation['POS'])
                bamsurgeon_snv_dict['vaf'].append(vaf.split(',')[mi])
                bamsurgeon_snv_dict['alt'].append(mut[0])
            else:  # mutation is insertion or deletion
                if mutation_type[mi] == 'DEL':
                    len_mut = len(mutation['REF']) - len(mut)  # important for BamSurgeon
                else:  # mutation_type[mi] == 'INS'
                    len_mut = len(mut) - len(mutation['REF'])  # not important for BamSurgeon
                bamsurgeon_indel_dict['chr'].append(str(mutation['#CHROM']))
                bamsurgeon_indel_dict['pos_start'].append(mutation['POS']-1)
                bamsurgeon_indel_dict['pos_end'].append(mutation['POS']+len_mut-1)
                bamsurgeon_indel_dict['vaf'].append(vaf.split(',')[mi])
                bamsurgeon_indel_dict['type'].append(mutation_type[mi])
                bamsurgeon_indel_dict['alt'].append(mut if mutation_type[mi] == 'INS' else '')
    return bamsurgeon_snv_dict, bamsurgeon_indel_dict


def get_mutation_type(ref, alt):
    if ',' not in alt:
        if len(ref) == len(alt):  # mutation_type == 'SNV'
            if len(alt) != 1:
                print('lengths of ref {} and alt {} are equal but not single nucleotide base'.format(ref, alt))
                if not (alt[0] != ref[0] and alt[1:] == ref[1:]):
                    raise ValueError('lengths of ref {} and alt {} are equal but not single nucleotide base change'.format(ref, alt))
            mutation_type = 'SNV'
        elif len(ref) < len(alt):
            mutation_type = 'INS'
        elif len(ref) > len(alt):
            mutation_type = 'DEL'
        else:
            raise ValueError('cannot find mutation type of REF={}, ALT={}'.format(ref, alt))
    else:  # several variants
        mutation_type = []
        for mi, mut in enumerate(alt.split(',')):
            if len(ref) == len(mut):  # mutation_type == 'SNV'
                if len(mut) != 1:
                    print('lengths of ref {} and alt {} are equal but not single nucleotide base'.format(ref, mut))
                    if not (mut[0] != ref[0] and mut[1:] == ref[1:]):
                        raise ValueError('lengths of ref {} and alt {} are equal but not single nucleotide base change'.format(ref, mut))
                mutation_type.append('SNV')
            elif len(ref) < len(mut):
                mutation_type.append('INS')
            elif len(ref) > len(mut):
                mutation_type.append('DEL')
            else:
                raise ValueError('cannot find mutation type of REF={}, ALT={}'.format(ref, mut))
    return mutation_type


def get_listpatientsnps(patient_snps_df):
    # list positions in patient
    old_listpatientssnps = list(patient_snps_df['ID'])
    old_listpatientssnps = [i for i in old_listpatientssnps if i != '.'] # remove unknown snps
    # multiple rsids at the same locus => replace 'rsid1;rsid2' by 'rsid1', 'rsid2'
    listpatientssnps = []
    for ki in old_listpatientssnps:
        if ';' in ki:
            for li in ki.split(';'):
                listpatientssnps.append(li)
        elif ki != '.':
            listpatientssnps.append(ki)
    return listpatientssnps




if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import seaborn as sns

    bamfile_path = '../data/healthy_chr22_merged-ready.bam'
    vcf_df = pd.read_csv('../data/common_SNPs/dbsnp_df.csv')
    patient_snps = pd.read_csv('../data/patient_SNPs/patient_snps.csv')

    reads2remove, log_pd = list_reads_to_remove(bamfile_path, vcf_df.iloc[3700:3750], patient_snps, verbose=1)
    '''
    plt.figure(figsize=(10, 5))
    plt.title('SNV')
    sns.histplot(data=log_pd[log_pd['type'] == 'SNV'][['vaf', 'normal af', 'noisy af']], bins=100,  stat="probability")
    plt.figure(figsize=(10, 5))
    plt.title('DEL')
    sns.histplot(data=log_pd[log_pd['type'] == 'DEL'][['vaf', 'normal af', 'noisy af']], bins=100,  stat="probability")
    plt.figure(figsize=(10, 5))
    plt.title('INS')
    sns.histplot(data=log_pd[log_pd['type'] == 'INS'][['vaf', 'normal af', 'noisy af']], bins=100,  stat="probability")
    plt.show()
    '''

    bamsurgeon_snv_pd, bamsurgeon_indel_pd = prepare_bamsurgeon_inputs(patient_snps, log_pd, max_vaf=0.1)
    print(bamsurgeon_snv_pd.head(10))
    print(bamsurgeon_indel_pd.head(10))


