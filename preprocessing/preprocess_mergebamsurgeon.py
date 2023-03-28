import os
import pysam
import argparse
from uuid import uuid4
from shutil import move
import logging

from bamsurgeon.common import *
from bamsurgeon.markreads import markreads
import bamsurgeon.replacereads as rr

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--bamfilename', dest='bamFileName', required=True, help='initial bam file')
parser.add_argument('-t', '--tmpbamfolder', dest='tmp_bam_folder', required=True, help='sam/bam file from which to obtain reads')
parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
parser.add_argument('-m', '--mutsfolder', dest='outbam_mutsfolder', required=True, help='folder containing all the temporaty mutation files to be merged')
parser.add_argument('--tagreads', action='store_true', default=False, help='add BS tag to altered reads')
parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
parser.add_argument('--seed', default=None, help='seed random number generation')

args = parser.parse_args()
print(args)

FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def replace(origbamfile, mutbamfile, outbamfile, seed=None):
    ''' open .bam file and call replacereads
    '''
    origbam = pysam.Samfile(origbamfile, 'rb')
    mutbam = pysam.Samfile(mutbamfile, 'rb')
    outbam = pysam.Samfile(outbamfile, 'wb', template=origbam)

    rr.replaceReads(origbam, mutbam, outbam, keepqual=True, seed=seed)

    origbam.close()
    mutbam.close()
    outbam.close()


# make a temporary file to hold mutated reads
outbam_mutsfile = "addsnv." + str(uuid4()) + ".muts.bam"
bamfile = pysam.Samfile(args.bamFileName, 'rb')
outbam_muts = pysam.Samfile(outbam_mutsfile, 'wb', template=bamfile)
outbam_muts.close()
bamfile.close()
print(outbam_mutsfile)

tmpbams = [os.listdir(tpb) for tpb in args.tmp_bam_folder if tpb.endswith('.bam') and tpb.startswith('add')]
if len(tmpbams) <= 1:
    raise ValueError('only {} successful mutations'.format(len(tmpbams)))
else:
    print('{} successful mutations'.format(print(len(tmpbams))))

# merge tmp bams
mergebams(tmpbams, outbam_mutsfile, maxopen=int(args.maxopen))

# cleanup
# for bam in tmpbams:
#    if os.path.exists(bam):
#        os.remove(bam)
#    if os.path.exists(bam + '.bai'):
#        os.remove(bam + '.bai')

if args.tagreads:
    tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
    print(tmp_tag_bam)
    markreads(outbam_mutsfile, tmp_tag_bam)
    move(tmp_tag_bam, outbam_mutsfile)
    logger.info("tagged reads.")

logger.info("done making mutations, merging mutations into %s --> %s" % (args.bamFileName, args.outBamFile))
replace(args.bamFileName, outbam_mutsfile, args.outBamFile, seed=args.seed)

# cleanup
# os.remove(args.outbam_mutsfile)

