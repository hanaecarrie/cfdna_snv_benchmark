#!/usr/bin/python

import argparse
import os
import json
import subprocess
import sys


# Example:

# python generate_fastqs.py -i <path to bam file> -t <tmp dir, will be created if not exists> -o <output directory, will be created if not exists>


# python generate_fastqs.py -i /mnt/projects/carriehc/cfDNA/cfdna_snv/cfdna_snv_benchmark/data/dilutions_series_chr22/dilutions_CRC-986_100215/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0/dilution_chr22_CRC-986_100215_1_CRC-986_300316_0.sorted.bam -t /mnt/projects/skanderupamj/wgs/data/raw.fastqs/test-picard-2-fq/tmp -o /mnt/projects/skanderupamj/wgs/data/raw.fastqs/test-picard-2-fq/results


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--inbam', help='Input bam file from which fastqs will be generated')
    parser.add_argument('-t', '--tmpdir', help='Temporary work directory, will be created if does not exist.')
    parser.add_argument('-o', '--outdir', help='Output directory where fastqs will be placed')
    parser.add_argument('-d', '--debug', default=False, help='Debug mode, tmp dir will not be deleted')
    parser.add_argument('-c', '--config', default="fastq.extraction.json", help='Config file showing paths to tools')

    args = parser.parse_args()

    is_debug = args.debug
    inbam = args.inbam
    tmpdir = args.tmpdir
    try:
        os.stat(tmpdir)
    except:
        os.makedirs(tmpdir)
        
    outdir = args.outdir

    # We assume the bam file ends with bam
    filename = os.path.basename(inbam)[:-4]

    unfiltered_sam = os.path.join(tmpdir, "%s.unfiltered.sam" %(filename))
    filtered_sam = os.path.join(tmpdir, "%s.filtered.sam" %(filename))
    filtered_bam = os.path.join(tmpdir, "%s.filtered.bam" %(filename))
    fq1 = os.path.join(outdir, "%s_R1.fastq.gz" %(filename))
    fq2 = os.path.join(outdir, "%s_R2.fastq.gz" %(filename))    
    
    configfile = args.config
    config = {}
    with open(configfile, 'r') as cfile:
        config = json.load(cfile)

    print(config)

    samtools = None
    if "samtools" in config:
        samtools = config["samtools"]
    java = None
    if "java" in config:
        java = config["java"]

    picard_jar = None
    if "picard" in config:
        picard_jar = config["picard"]
    

    if is_debug is False and outdir == tmpdir:
        print("Do not specify output bam file in tmp directory, as the tmp directory will be deleted by this script")
        sys.exit(1)


    # Step 1: convert bam file to sam
    print("Generating SAM file")
    cmd = "samtools view -h %s > %s" %(inbam, unfiltered_sam)
    if samtools is not None:
        cmd = "%s view -h %s > %s" %(samtools, inbam, unfiltered_sam)
    try:
        print(cmd)
        subprocess.call(cmd, shell = True)
    except subprocess.CalledProcessError as e:
        print("BAM to SAM error")
        print(e)
        sys.exit(1)

    # Step 2: filter out unpaired reads from sam
    reads = {} # K00115:204:H7G7GBBXX:1:1111:14387:45379 -> cnt


    print("Identifying unpaired reads from SAM file")
    with open(unfiltered_sam, 'r') as ifile:
        for row in ifile:
            if not row.startswith("@"):
                dat = row.rstrip('\n').strip().split('\t')
                identity = dat[0]
                if not identity in reads:
                    reads[identity] = 0

                cnt = reads[identity]
                cnt += 1
                reads[identity] = cnt

    paired_reads = {}
    unpaired_reads = {}
    strange_reads = {}

    for identity in reads:
        cnt = reads[identity]
        if cnt == 2:
            paired_reads[identity] = 0
        elif cnt == 1:
            unpaired_reads[identity] = 0
        else:
            strange_reads[identity] = 0


    paired = float(len(paired_reads))
    unpaired = float(len(unpaired_reads))
    strange = float(len(strange_reads))
    total = paired + unpaired + strange

    pct_unpaired = unpaired / total * 100.0
    
            
    print("Number of paired reads: %d" %(len(paired_reads)))
    print("Number of unpaired reads: %d" %(len(unpaired_reads)))
    print("Number of reads with more than 2 reads: %d" %(len(strange_reads)))
    print("Percentage of unpaired reads = %s percent" %(str(pct_unpaired)))
    print("\n\n")
    
    print("Creating SAM file with only paired reads")
    with open(unfiltered_sam, 'r') as ifile, open(filtered_sam, 'w') as ofile:
        for row in ifile:
            if row.startswith("@"):
                ofile.write(row)
            else:
                dat = row.rstrip('\n').strip().split('\t')
                identity = dat[0]
                if identity in paired_reads:
                    ofile.write(row)

    # Step 3: Generate bam from sam
    print("Creating BAM file from SAM file")
    cmd = "samtools view -b %s > %s" %(filtered_sam, filtered_bam)
    if samtools is not None:
        cmd = "%s view -b %s > %s" %(samtools, filtered_sam, filtered_bam)

    
    try:
        print(cmd)
        subprocess.call(cmd, shell = True)
    except subprocess.CalledProcessError as e:
        print("Error when generating bam file from filtered sam file")
        print(e)
        sys.exit(1)

    # Step 4: Extracting paired fastqs from BAM file
    print("Extracting fastqs")
    print("Making output directory...if not exists.")
    try:
        os.stat(outdir)
    except:
        os.makedirs(outdir)

        
    cmd = "%s -Xmx2g -jar %s SamToFastq -I %s -F %s -F2 %s" %(java, picard_jar, filtered_bam, fq1, fq2)
    try:
        print(cmd)
        subprocess.call(cmd, shell = True)
    except subprocess.CalledProcessError as e:
        print("Error extracting fastqs")
        print(e)
        sys.exit(1)


    # Step 5: Cleaning up
    if is_debug is False:
        cmd = "rm -rf %s" %(tmpdir)
        try:
            subprocess.call(cmd, shell = True)
        except subprocess.CalledProcessError as e:
            print("Unable to delete tmp directory: %s" %(tmpdir))
            print(e)
            sys.exit(1)


        
    print("Completed:")
    print(fq1)
    print(fq2)
    
        
    
        
    


#------------------------------------------------------

if __name__ == "__main__":
    main()


    
