#! /usr/bin/env python
import argparse
import fastq_files as ff


parser = argparse.ArgumentParser(description='add required annotations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dir',type=str,help='Reference file',required=True)
parser.add_argument('--r1',type=str,help='Sample name',required=True)
parser.add_argument('--r2',type=str,help='Sample name',required=True)

args = parser.parse_args()


for sample in ff.find_fastq_files([args.dir],args.r1,args.r2):
    print(sample)