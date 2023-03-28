#!/usr/bin/python

import numpy as np
from pyliftover import LiftOver


def liftover(coords, chainfile):
    # hg19to38: chr7_1920403
    # coords = chrom_startpos_endpos or chrom_pos
    lo = LiftOver(chainfile)
    sc = coords.split('_')
    if len(sc) == 3:
        chrom, startpos, endpos = coords.split('_')
    else:  # len(sc) == 2:
        chrom, startpos = coords.split('_')
    startpos = int(startpos)
    startpos_new = lo.convert_coordinate(chrom, startpos)
    if len(sc) == 3:
        endpos = int(endpos)
        endpos_new = lo.convert_coordinate(chrom, endpos)
        if (startpos_new != np.nan) and (endpos_new != np.nan) and (startpos_new != []) and (endpos_new != []):
            chrom_new = startpos_new[0][0].replace("chr", "")
            startpos_new = startpos_new[0][1]
            endpos_new = endpos_new[0][1]
            coords_new = str(chrom_new) + '_' + str(startpos_new) + '_' + str(endpos_new)
            return coords_new
        else:
            return 'nan_nan_nan'
    else:  # len(sc) == 2:
        if (startpos_new != np.nan) and (startpos_new != []):
            chrom_new = startpos_new[0][0].replace("chr", "")
            startpos_new = startpos_new[0][1]
            coords_new = str(chrom_new) + '_' + str(startpos_new)
            return coords_new
        else:
            print('fail')
            return 'nan_nan'
