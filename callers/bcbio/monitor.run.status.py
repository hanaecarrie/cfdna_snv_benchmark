#!/usr/bin/python

import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', help="Directory of bcbio runs")
    args = parser.parse_args()

    rootdir = args.directory
    samples = os.listdir(rootdir)

    for sample in samples:
        logfile = os.path.join(rootdir, sample, "bcbio_work", "z.%s.e" %(sample))
        if not os.path.exists(logfile):
            logfile = os.path.join(rootdir, sample, "work", "log", "bcbio-nextgen.log")

        if os.path.exists(logfile):
            try:
                process = subprocess.Popen(["tail", "-n1", logfile], stdout=subprocess.PIPE)
                result = process.communicate()[0]
                result = result.decode('utf-8')
                result = result.rstrip('\n').strip()
                print("%s\t%s" % (sample, result))
            except subprocess.CalledProcessError as e:
                print(e)


if __name__ == "__main__":
    main()
    

