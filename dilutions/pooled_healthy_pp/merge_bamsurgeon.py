
import sys

if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(outbam_mutsfile, tmp_tag_bam)
            move(tmp_tag_bam, outbam_mutsfile)
            logger.info("tagged reads.")

        logger.info("done making mutations, merging mutations into %s --> %s" % (args.bamFileName, args.outBamFile))
        replace(args.bamFileName, outbam_mutsfile, args.outBamFile, seed=args.seed)

        #cleanup
        os.remove(outbam_mutsfile)
