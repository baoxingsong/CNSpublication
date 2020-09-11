#!python
import re
import subprocess
import sys
from argparse import ArgumentParser
import sys
#read a fasta file and return a dictionary, the key is entry id and the value is the sequence in upcase
from utils import readFastaFile
from utils import str2bool
import re

class SNP:
    def __init__(self, chr, v4cordinate, maf, geno, depth, mpileup):
        self.chr = chr
        self.v4cordinate = v4cordinate
        self.maf = maf
        self.geno = geno
        self.depth = depth
        self.mpileup = mpileup
    def __str__(self):
        return (self.chr + "\t" + str(self.v4cordinate) + "\t" + str(self.maf) + "\t" + str(self.geno) + "\t" + str(self.depth) + "\t" + str(self.mpileup))


if __name__ == '__main__':
    parser = ArgumentParser(description='count number of based overlap between CNS bam output and the eQTL result,'
                                        'please input the vcf file and the eqtl for one chromosome only')

    parser.add_argument("-g", "--genome",
                        dest="genome",
                        type=str,
                        default="",
                        help="the masked reference genome file")

    parser.add_argument("-b", "--bam",
                        dest="bam",
                        type=str,
                        default="",
                        help="the output of and-CNS pipeline in bam format")

    parser.add_argument("-c", "--chr",
                        dest="chr",
                        type=str,
                        default="",
                        help="the chromosome to be analysised")

    parser.add_argument("-s", "--mask",
                        dest="mask",
                        type=str2bool,
                        default=True,
                        help="only count the non-masking region SNP and genome length")

    parser.add_argument("-m", "--miss",
                        dest="miss",
                        type=str,
                        default="",
                        help="the missing statistics using V4 coordinate")

    parser.add_argument("-v", "--vcf",
                        dest="vcf",
                        type=str,
                        default="",
                        help="the B73 v4 variant file in vcf format")

    parser.add_argument("-f", "--maf",
                        dest="maf",
                        type=str,
                        default="",
                        help="the MAF statistics using V4 coordinate")

    args = parser.parse_args()


    if args.genome == "":
        print("Error: please specify --genome", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.bam == "":
        print("Error: please specify --bam", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.miss == "":
        print("Error: please specify --miss", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.maf == "":
        print("Error: please specify --maf", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.vcf == "":
        print("Error: please specify --vcf", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.chr == "":
        print("Error: please specify --chr", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    reference_genome = readFastaFile(args.genome)
    print("reference genome reading done", file=sys.stderr)

    chr = args.chr
    id_to_v4cordinate_dict = dict()
    with open(args.vcf) as f:
        for line in f:
            if line[0] is not '#':
                elements = line.split('\t')
                if elements[0] == chr:
                    id_to_v4cordinate_dict[elements[2]] = elements[1]


    sig_v4cordinate_dict = dict()
    totalDepth = 0
    totalMpileup = 0
    # print("SNP reading done", file=sys.stderr)
    seq = reference_genome[chr]
    seq = re.sub("\\s", "", seq)
    seq = re.sub("-", "", seq)
    if args.mask:
        # read the VCF file begin
        with open(args.maf) as f:
            for line in f:
                line = line.lstrip()
                if line[0] is not 'C':
                    elements = line.split()
                    position = id_to_v4cordinate_dict[elements[1]]
                    if chr == elements[0] and (reference_genome[chr][int(position)-1] is not 'n') and (elements[4][0] is not 'N'):
                        s = SNP(elements[0], int(position), float(elements[4]), -1, 0, 0)
                        sig_v4cordinate_dict[position] = s
        # print("vcf file reading done", file=sys.stderr)
        with open(args.miss) as f:
            for line in f:
                line = line.lstrip()
                if line[0] is not 'C':
                    elements = line.split()
                    position = id_to_v4cordinate_dict[elements[1]]
                    if chr == elements[0] and (reference_genome[chr][int(position)-1] is not 'n') and (position in sig_v4cordinate_dict) and (elements[4][0] is not 'N'):
                        sig_v4cordinate_dict[position].geno = float(elements[4])

        for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split()
                position = elements2[1]
                if int(elements2[2]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n'):
                    totalDepth = totalDepth + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].depth = int(elements2[2])
                    # print()
                    # print(line2)
                    # print(sig_v4cordinate_dict[position])
        # print("samtools depth done", file=sys.stderr)

        for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split()
                position = elements2[1]
                if int(elements2[3]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n'):
                    totalMpileup = totalMpileup + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].mpileup = int(elements2[3])
        seq = re.sub("n", "", seq)
    else:
        with open(args.maf) as f:
            for line in f:
                line = line.lstrip()
                if line[0] is not 'C':
                    elements = line.split()
                    position = id_to_v4cordinate_dict[elements[1]]
                    if (chr == elements[0]  and (elements[4][0] is not 'N')):
                        s = SNP(elements[0], int(position), float(elements[4]), -1, 0, 0)
                        sig_v4cordinate_dict[position] = s
            # print("vcf file reading done", file=sys.stderr)
        with open(args.miss) as f:
            for line in f:
                line = line.lstrip()
                if line[0] is not 'C':
                    elements = line.split()
                    position = id_to_v4cordinate_dict[elements[1]]
                    if (chr == elements[0]  and (elements[4][0] is not 'N') and (position in sig_v4cordinate_dict) ):
                        sig_v4cordinate_dict[position].geno = float(elements[4])

        for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split()
                position = elements2[1]
                if int(elements2[2]) > 0:
                    totalDepth = totalDepth + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].depth = int(elements2[2])
        # print("samtools depth done", file=sys.stderr)

        for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split()
                position = elements2[1]
                if int(elements2[3]) > 0:
                    totalMpileup = totalMpileup + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].mpileup = int(elements2[3])

    for position in sig_v4cordinate_dict:
        if sig_v4cordinate_dict[position].mpileup > 0:
            print("mpileup\t" + str(sig_v4cordinate_dict[position].geno) + "\t" + str(sig_v4cordinate_dict[position].maf))
        else:
            print("non-mpileup\t" + str(sig_v4cordinate_dict[position].geno) + "\t" + str(sig_v4cordinate_dict[position].maf))
        if sig_v4cordinate_dict[position].depth > 0:
            print("depth\t" + str(sig_v4cordinate_dict[position].geno) + "\t" + str(sig_v4cordinate_dict[position].maf))
        else:
            print("non-depth\t" + str(sig_v4cordinate_dict[position].geno) + "\t" + str(sig_v4cordinate_dict[position].maf))
        # print(sig_v4cordinate_dict[position])
