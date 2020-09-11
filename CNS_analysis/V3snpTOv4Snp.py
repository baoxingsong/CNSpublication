#!python
import re
import subprocess
import sys
from argparse import ArgumentParser
import sys
#read a fasta file and return a dictionary, the key is entry id and the value is the sequence in upcase
from utils import readFastaFile
from utils import SNP
from utils import str2bool
import re

if __name__ == '__main__':
    parser = ArgumentParser(description='transform the VS genotypic variants from V3 to V4')


    parser.add_argument("-m", "--hapmap",
                        dest="hapmap",
                        type=str,
                        default="",
                        help="hapmap file")

    parser.add_argument("-v", "--vcf",
                        dest="vcf",
                        type=str,
                        default="",
                        help="the B73 v4 variant file in vcf format")
    args = parser.parse_args()


    if args.hapmap == "":
        print("Error: please specify --hapmap", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.vcf == "":
        print("Error: please specify --vcf", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    chr = ""
    snps_dict = dict()
    sig_v4cordinate_dict = dict()

    # read the SNP hapmap begin
    with open(args.hapmap) as f:
        for line in f:
            if line[0] is 'S':
                elements = line[:100].split('\t')
                chr = elements[2]
                s = SNP(elements[2], int(elements[3]), int(elements[3]), 0, 0)
                snps_dict[elements[0]] = s
    # read the SNP hapmap end
    # print("SNP reading done", file=sys.stderr)
    with open(args.vcf) as f:
        for line in f:
            if line[0] is not '#':
                elements = line.split('\t')
                if elements[0] == chr:
                    elements2 = elements[2].split('-')
                    variant_id = "S" + elements2[0] + "_" + elements2[1]
                    if (variant_id in snps_dict):
                        print(chr + "\t" + str(elements[1]))
    # print("vcf file reading done", file=sys.stderr)
