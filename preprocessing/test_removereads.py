import unittest
import random
import pysam
import re
import pandas as pd

from supporting_reads import get_mutation_type, assess_mutation


class MyTestCase(unittest.TestCase):
    nucleotides_ref = ['A', 'T', 'C', 'G']
    nucleotides_alt = ['A', 'T', 'C', 'G', 'N']

    def test_mutation_type(self):

        print('### test input values ###')
        ref, alt, = 'NTG', 'TG'
        self.assertRaises(ValueError, get_mutation_type, ref, alt)
        ref, alt, = 'TG', 'TGa'
        self.assertRaises(ValueError, get_mutation_type, ref, alt)
        ref, alt, = 'TG', 'TG,TGB'
        self.assertRaises(ValueError, get_mutation_type, ref, alt)

        print('### test reference ###')
        # reference
        ref = random.choice(self.nucleotides_ref)
        alt = ref
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'ref'))
        self.assertEqual(get_mutation_type(ref, alt), 'ref')

        print('### test substitution ###')
        # substitution, SNV with 1 nucleotide
        ref = random.choice(self.nucleotides_ref)
        alt = random.choice([i for i in self.nucleotides_alt if i != ref])
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'SNV'))
        self.assertEqual(get_mutation_type(ref, alt), 'SNV')

        # substitution, SNV with several nucleotides
        lenmut = random.randint(1, 20)
        initialnucleotide = random.choice(self.nucleotides_ref)
        ref = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = random.choice([i for i in self.nucleotides_ref if i != initialnucleotide]) + ref[1:]
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'SNV'))
        self.assertEqual(get_mutation_type(ref, alt), 'SNV')

        # other, looks like SNV with several nucleotides but other changes in sequence
        lenmut = random.randint(1, 20)
        initialnucleotide = random.choice(self.nucleotides_ref)
        ref = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = random.choice([i for i in self.nucleotides_ref if i != initialnucleotide]) + ref[1:-1] + 'N'
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'other'))
        self.assertEqual(get_mutation_type(ref, alt), 'other')

        print('### test deletion ###')
        # deletion
        lenmut = random.randint(2, 20)
        initialnucleotide = random.choice(self.nucleotides_ref)
        ref = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = initialnucleotide
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'DEL'))
        self.assertEqual(get_mutation_type(ref, alt), 'DEL')

        ref = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = initialnucleotide + random.choice(self.nucleotides_alt)
        self.assertEqual(get_mutation_type(ref, alt), 'DEL')

        # unknown, looks like deletion but different initial nucleotide
        ref = random.choice([i for i in self.nucleotides_ref if i != initialnucleotide]) + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = initialnucleotide
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'other'))
        self.assertEqual(get_mutation_type(ref, alt), 'other')

        print('### test insertion ###')
        # insertion
        lenmut = random.randint(2, 20)
        initialnucleotide = random.choice(self.nucleotides_ref)
        ref = initialnucleotide
        alt = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'INS'))
        self.assertEqual(get_mutation_type(ref, alt), 'INS')

        ref = initialnucleotide + random.choice(self.nucleotides_ref)
        alt = initialnucleotide + ''.join([random.choice(self.nucleotides_alt) for _ in range(lenmut)])
        self.assertEqual(get_mutation_type(ref, alt), 'INS')

        # unknown, looks like deletion but different initial nucleotide
        ref = initialnucleotide
        alt = random.choice([i for i in self.nucleotides_alt if i != initialnucleotide]) + ''.join([random.choice(self.nucleotides_alt) for _ in range(lenmut)])
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'other'))
        self.assertEqual(get_mutation_type(ref, alt), 'other')

        print('### test list ###')
        # list
        ref = 'T'
        alt = 'TCA,A,AG'
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), ['INS', 'SNV', 'other']))
        self.assertEqual(get_mutation_type(ref, alt), ['INS', 'SNV', 'other'])

    def test_assess_mutation(self):

        vcf_df = pd.read_csv('../data/common_SNPs/dbsnp_df.csv')

        snv_done, ins_done, del_done = False, False, False

        while not (snv_done and ins_done and del_done):
            print(snv_done + ins_done + del_done)

            rand_mut = random.randint(0, vcf_df.shape[0])
            mutation = vcf_df.iloc[rand_mut]

            genotype = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}  # initiialise genotype (for SNVs only)

            bamfile_path = '../data/healthy_chr22_merged-ready.bam'
            samfile = pysam.AlignmentFile(bamfile_path, "rb")
            print(str(mutation['#CHROM']), mutation['POS']-1, mutation['POS'])
            for pileupcolumn in samfile.pileup(str(mutation['#CHROM']), mutation['POS']-1, mutation['POS'], min_base_quality=0):
                if pileupcolumn.pos == mutation['POS']-1:
                    print("\ndepth of coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                    depthcov = pileupcolumn.n
            rand_read = random.randint(0, depthcov+1)  # cov, randint half open intervals [a,b)
            index_reads = [0, 1, rand_read, depthcov-1, depthcov]
            print(rand_read)
            count_read = 0

            for read in samfile.fetch(str(mutation['#CHROM']), mutation['POS']-1, mutation['POS']):
                if count_read in index_reads:
                    cond, mutation_type, seq, pos, cigar, genotype = assess_mutation(read, mutation, genotype)
                    if mutation_type == 'SNV':
                        print(cond, mutation_type, mutation['REF'], mutation['ALT'], seq[pos], pos, cigar, genotype)
                        if cond:
                            snv_done = True
                            self.assertEqual(mutation['ALT'], seq[pos])
                        else:
                            self.assertNotEqual(mutation['ALT'], seq[pos])
                    else:
                        print(cond, mutation_type, mutation['REF'], mutation['ALT'], seq[pos:pos+max(len(mutation['REF']), len(mutation['ALT']))], pos, cigar)
                        if mutation_type == 'INS' and cond:
                            ins_done = True
                            cigar_states = re.split('[0-9]+', cigar)[1:]
                            cigar_pos = re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
                            cigar_pos = [0 if (cigar_states[i] == 'S') else int(ci) for i, ci in enumerate(cigar_pos)]
                            self.assertEqual(sum(cigar_pos[:cigar_states.index('I')]), pos+1)
                        elif mutation_type == 'DEL' and cond:
                            del_done = True
                            cigar_states = re.split('[0-9]+', cigar)[1:]
                            cigar_pos = re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
                            cigar_pos = [0 if (cigar_states[i] == 'S') else int(ci) for i, ci in enumerate(cigar_pos)]
                            self.assertEqual(sum(cigar_pos[:cigar_states.index('D')]), pos+1)
                else:
                    pass
                count_read += 1


if __name__ == '__main__':
    unittest.main()
