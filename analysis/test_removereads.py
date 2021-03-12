import unittest
import random

from supporting_reads import get_mutation_type


class MyTestCase(unittest.TestCase):
    nucleotides_ref = ['A', 'T', 'C', 'G']
    nucleotides_alt = ['A', 'T', 'C', 'G', 'N']

    def test_mutation_type(self):

        ref, alt, = 'NTG', 'TG'
        self.assertRaises(ValueError, get_mutation_type, ref, alt)
        ref, alt, = 'TG', 'TGa'
        self.assertRaises(ValueError, get_mutation_type, ref, alt)
        ref, alt, = 'TG', 'TG,TGB'
        self.assertRaises(ValueError, get_mutation_type, ref, alt)

        # reference
        ref = random.choice(self.nucleotides_ref)
        alt = ref
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'ref'))
        self.assertEqual(get_mutation_type(ref, alt), 'ref')

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

        ###
        #print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'other'))
        #self.assertEqual(get_mutation_type(ref, alt), 'other')

        # deletion
        lenmut = random.randint(1, 20)
        initialnucleotide = random.choice(self.nucleotides_ref)
        ref = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = initialnucleotide
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'DEL'))
        self.assertEqual(get_mutation_type(ref, alt), 'DEL')

        ref = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = initialnucleotide + random.choice(self.nucleotides_alt)
        self.assertRaises(ValueError, get_mutation_type, ref, alt)

        ref = random.choice([i for i in self.nucleotides_ref if i != initialnucleotide]) + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        alt = initialnucleotide
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'unknown'))
        self.assertEqual(get_mutation_type(ref, alt), 'unknown')

        # insertion
        lenmut = random.randint(1, 20)
        initialnucleotide = random.choice(self.nucleotides_ref)
        ref = initialnucleotide
        alt = initialnucleotide + ''.join([random.choice(self.nucleotides_ref) for _ in range(lenmut)])
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'INS'))
        self.assertEqual(get_mutation_type(ref, alt), 'INS')

        ref = initialnucleotide + random.choice(self.nucleotides_ref)
        alt = initialnucleotide + ''.join([random.choice(self.nucleotides_alt) for _ in range(lenmut)])
        self.assertRaises(ValueError, get_mutation_type, ref, alt)

        ref = initialnucleotide
        alt = random.choice([i for i in self.nucleotides_alt if i != initialnucleotide]) + ''.join([random.choice(self.nucleotides_alt) for _ in range(lenmut)])
        print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'unknown'))
        self.assertEqual(get_mutation_type(ref, alt), 'unknown')

        ###
        #print('REF: {}, ALT: {}, detected mutation type: {}, expected mutation type: {}'.format(ref, alt, get_mutation_type(ref, alt), 'unknown'))
        #self.assertEqual(get_mutation_type(ref, alt), 'unknown')

    #def test_mutation_type_list(self):
    #    lenmuts_ref = [random.randint(1, 5) for _ in range(4)]
    #    ref = [random.choice(self.nucleotides) for _ in range(lenmut) for lenmut in lenmuts_ref]
    #    alt = [random.choice(self.nucleotides) for _ in range(lenmut) for lenmut in lenmuts_alt]


if __name__ == '__main__':
    unittest.main()
