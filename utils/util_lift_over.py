#!/usr/bin/python

import os
from pyliftover import LiftOver

def liftMeOver(chainfile, infile, outfile):


    lo = LiftOver(chainfile)
    with open(infile, 'r') as ifile, open(outfile, 'w') as ofile:
        ofile.write(ifile.readline())

        for row in ifile:
            row = row.rstrip('\n').strip()
            data = row.split('\t')

            chrom = data[1]
            spos = int(data[2])
            epos = int(data[3])

            spos_hg38 = lo.convert_coordinate(chrom, spos)
            epos_hg38 = lo.convert_coordinate(chrom, epos)
    
            if len(spos_hg38) > 0 and len(epos_hg38):
                start_tuple = spos_hg38[0]
                end_tuple = epos_hg38[0]

                sposition = start_tuple[1]
                eposition = end_tuple[1]

                if sposition < eposition:
                
                    data[2] = str(sposition)
                    data[3] = str(eposition)

                else:
                    data[2] = str(eposition)
                    data[3] = str(sposition)

                    
                # I need to remove chr
                chrom = chrom.replace("chr", "")
                
                nrow = "\t".join(data)
                ofile.write("%s\n" %(nrow))


    
