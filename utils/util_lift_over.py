#!/usr/bin/python

import numpy as np
from pyliftover import LiftOver


def liftover(coords, chainfile):
    # coords = chrom_startpos_endpos
    lo = LiftOver(chainfile)
    chrom, startpos, endpos = coords.split('_')
    startpos = int(startpos)
    endpos = int(endpos)
    startpos_new = lo.convert_coordinate(chrom, startpos)
    endpos_new = lo.convert_coordinate(chrom, endpos)
    if (startpos_new != np.nan) and (endpos_new != np.nan) and (startpos_new != []) and (endpos_new != []):
        chrom_new = startpos_new[0][0].replace("chr", "")
        startpos_new = startpos_new[0][1]
        endpos_new = endpos_new[0][1]
        coords_new = str(chrom_new) + '_' + str(startpos_new) + '_' + str(endpos_new)
        return coords_new
    else:
        return 'nan_nan_nan'
