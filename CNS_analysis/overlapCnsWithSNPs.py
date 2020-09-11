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
    def __init__(self, chr, v4cordinate, depth, mpileup):
        self.chr = chr
        self.v4cordinate = v4cordinate
        self.depth = depth
        self.mpileup = mpileup
    def __str__(self):
        return (self.chr + "\t" + "\t" + str(self.v4cordinate) + "\t" + str(self.depth) + "\t" + str(self.mpileup))


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

    parser.add_argument("-v", "--vcf",
                        dest="vcf",
                        type=str,
                        default="",
                        help="the B73 v4 variant file in vcf format")

    parser.add_argument("-m", "--bim",
                        dest="bim",
                        type=str,
                        default="",
                        help="the B73 v4 variant file in plink bim format")
    args = parser.parse_args()


    if args.genome == "":
        print("Error: please specify --genome", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.bam == "":
        print("Error: please specify --bam", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.vcf == "" and args.bim == "":
        print("Error: please specify --vcf or --bim", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.chr == "":
        print("Error: please specify --chr", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    reference_genome = readFastaFile(args.genome)
    print("reference genome reading done", file=sys.stderr)

    totalDepth = 0
    totalMpileup = 0

    chr = args.chr
    snps = dict()

    # print("SNP reading done", file=sys.stderr)
    seq = reference_genome[chr]
    seq = re.sub("\\s", "", seq)
    seq = re.sub("-", "", seq)
    print ("chr" + chr)
    genomelength = len(reference_genome[chr])
    if args.mask:
        seq = seq.replace("n", "")
        seq = seq.replace("N", "")
        seq = seq.replace("b", "")
        seq = seq.replace("B", "")
        # read the VCF file begin
        if args.vcf != "" :
            with open(args.vcf) as f:
                for line in f:
                    if line[0] is not '#':
                        elements = line.split('\t')
                        if  (chr == elements[0] and (reference_genome[chr][int(elements[1])-1] is not 'n') and (reference_genome[chr][int(elements[1])-1] is not 'b') and (reference_genome[chr][int(elements[1])-1] is not 'N')  ) :
                            s = SNP(elements[0], int(elements[1]), 0, 0)
                            snps[elements[1]] = s
        else:
            with open(args.bim) as f:
                for line in f:
                    if line[0] is not '#':
                        elements = line.split('\t')
                        if  (chr == elements[0] and (reference_genome[chr][int(elements[1])-1] is not 'n') and (reference_genome[chr][int(elements[1])-1] is not 'b') and (reference_genome[chr][int(elements[1])-1] is not 'N')  ) :
                            s = SNP(elements[0], int(elements[3]), 0, 0)
                            snps[elements[3]] = s
        # print("vcf file reading done", file=sys.stderr)

        for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[2]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n') and (reference_genome[chr][int(elements2[1])-1] is not 'N') and (reference_genome[chr][int(elements2[1])-1] is not 'b') and (reference_genome[chr][int(elements2[1])-1] is not 'B'):
                    totalDepth = totalDepth + 1
                if position in snps:
                    snps[position].depth = int(elements2[2])
        # print("samtools depth done", file=sys.stderr)

        for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[3]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n') and (reference_genome[chr][int(elements2[1])-1] is not 'N') and (reference_genome[chr][int(elements2[1])-1] is not 'b') and (reference_genome[chr][int(elements2[1])-1] is not 'B'):
                    totalMpileup = totalMpileup + 1
                if position in snps:
                    snps[position].mpileup = int(elements2[3])
    else:
        if args.vcf != "":
            with open(args.vcf) as f:
                for line in f:
                    if line[0] is not '#':
                        elements = line.split('\t')
                        if chr == elements[0]:
                            s = SNP(elements[0], int(elements[1]), 0, 0)
                            snps[elements[1]] = s
        else:
            with open(args.bim) as f:
                for line in f:
                    if line[0] is not '#':
                        elements = line.split('\t')
                        if (chr == elements[0]):
                            s = SNP(elements[0], int(elements[3]), 0, 0)
                            snps[elements[3]] = s
        # print("vcf file reading done", file=sys.stderr)

        for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[2]) > 0:
                    totalDepth = totalDepth + 1
                if position in snps:
                    snps[position].depth = int(elements2[2])
        # print("samtools depth done", file=sys.stderr)

        for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[3]) > 0:
                    totalMpileup = totalMpileup + 1
                if position in snps:
                    snps[position].mpileup = int(elements2[3])

    number_snp_depth = 0
    number_snp_mpileup = 0
    for position in snps:
        if snps[position].depth > 0:
            number_snp_depth = number_snp_depth + 1
        if snps[position].mpileup > 0:
            number_snp_mpileup = number_snp_mpileup + 1

    print ("number_SNPs\t" + str(len(snps)))
    print ("totalDepth\t" + str(totalDepth))
    print ("totalMpileup\t" + str(totalMpileup))
    print ("totalChrLength\t" + str(len(seq)))
    print ("number_snp_depth\t" + str(number_snp_depth))
    print ("number_snp_mpileup\t" + str(number_snp_mpileup))




#
# if __name__ == '__main__':
#     parser = ArgumentParser(description='count number of based overlap between CNS bam output and the eQTL result,'
#                                         'please input the vcf file and the eqtl for one chromosome only')
#     parser.add_argument("-g", "--genome",
#                         dest="genome",
#                         type=str,
#                         default="",
#                         help="the masked reference genome file")
#
#     parser.add_argument("-b", "--bam",
#                         dest="bam",
#                         type=str,
#                         default="",
#                         help="the output of and-CNS pipeline in bam format")
#
#     parser.add_argument("-m", "--hapmap",
#                         dest="hapmap",
#                         type=str,
#                         default="",
#                         help="hapmap file")
#
#     parser.add_argument("-s", "--mask",
#                         dest="mask",
#                         type=str2bool,
#                         default=True,
#                         help="only count the non-masking region SNP and genome length")
#
#     parser.add_argument("-v", "--vcf",
#                         dest="vcf",
#                         type=str,
#                         default="",
#                         help="the B73 v4 variant file in vcf format")
#     args = parser.parse_args()
#
#
#     if args.genome == "":
#         print("Error: please specify --genome", file=sys.stderr)
#         parser.print_help()
#         sys.exit(1)
#
#     if args.bam == "":
#         print("Error: please specify --bam", file=sys.stderr)
#         parser.print_help()
#         sys.exit(1)
#
#     if args.hapmap == "":
#         print("Error: please specify --hapmap", file=sys.stderr)
#         parser.print_help()
#         sys.exit(1)
#
#     if args.vcf == "":
#         print("Error: please specify --vcf", file=sys.stderr)
#         parser.print_help()
#         sys.exit(1)
#
#     reference_genome = readFastaFile(args.genome)
#     print("reference genome reading done", file=sys.stderr)
#
#     totalDepth = 0
#     totalMpileup = 0
#
#     chr = ""
#     snps_dict = dict()
#     sig_v4cordinate_dict = dict()
#
#     # read the SNP hapmap begin
#     with open(args.hapmap) as f:
#         for line in f:
#             if line[0] is 'S':
#                 elements = line[:100].split('\t')
#                 chr = elements[2]
#                 s = SNP(elements[2], int(elements[3]), int(elements[3]), 0, 0)
#                 snps_dict[elements[0]] = s
#     # read the SNP hapmap end
#     # print("SNP reading done", file=sys.stderr)
#     seq = reference_genome[chr]
#     seq = re.sub("\\s", "", seq)
#     seq = re.sub("-", "", seq)
#     print ("chr" + chr)
#     if args.mask:
#         # read the VCF file begin
#         with open(args.vcf) as f:
#             for line in f:
#                 if line[0] is not '#':
#                     elements = line.split('\t')
#                     if elements[0] == chr:
#                         elements2 = elements[2].split('-')
#                         variant_id = "S" + elements2[0] + "_" + elements2[1]
#                         if (variant_id in snps_dict) and (reference_genome[chr][int(elements[1])-1] is not 'n'):
#                             sig_v4cordinate_dict[elements[1]] = snps_dict[variant_id]
#                             sig_v4cordinate_dict[elements[1]].v4cordinate = int(elements[1])
#         # print("vcf file reading done", file=sys.stderr)
#
#         for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
#             if len(line2) > 0:
#                 elements2 = line2.split('\t')
#                 position = elements2[1]
#                 if int(elements2[2]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n'):
#                     totalDepth = totalDepth + 1
#                 if position in sig_v4cordinate_dict:
#                     sig_v4cordinate_dict[position].depth = int(elements2[2])
#         # print("samtools depth done", file=sys.stderr)
#
#         for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
#             if len(line2) > 0:
#                 elements2 = line2.split('\t')
#                 position = elements2[1]
#                 if int(elements2[3]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n'):
#                     totalMpileup = totalMpileup + 1
#                 if position in sig_v4cordinate_dict:
#                     sig_v4cordinate_dict[position].mpileup = int(elements2[3])
#         seq = re.sub("n", "", seq)
#     else:
#         with open(args.vcf) as f:
#             for line in f:
#                 if line[0] is not '#':
#                     elements = line.split('\t')
#                     if elements[0] == chr:
#                         elements2 = elements[2].split('-')
#                         variant_id = "S" + elements2[0] + "_" + elements2[1]
#                         if (variant_id in snps_dict):
#                             sig_v4cordinate_dict[elements[1]] = snps_dict[variant_id]
#                             sig_v4cordinate_dict[elements[1]].v4cordinate = int(elements[1])
#         # print("vcf file reading done", file=sys.stderr)
#
#         for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
#             if len(line2) > 0:
#                 elements2 = line2.split('\t')
#                 position = elements2[1]
#                 if int(elements2[2]) > 0:
#                     totalDepth = totalDepth + 1
#                 if position in sig_v4cordinate_dict:
#                     sig_v4cordinate_dict[position].depth = int(elements2[2])
#         # print("samtools depth done", file=sys.stderr)
#
#         for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
#             if len(line2) > 0:
#                 elements2 = line2.split('\t')
#                 position = elements2[1]
#                 if int(elements2[3]) > 0:
#                     totalMpileup = totalMpileup + 1
#                 if position in sig_v4cordinate_dict:
#                     sig_v4cordinate_dict[position].mpileup = int(elements2[3])
#
#     number_eqtl_depth = 0
#     number_eqtl_mpileup = 0
#     for position in sig_v4cordinate_dict:
#         if sig_v4cordinate_dict[position].depth > 0:
#             number_eqtl_depth = number_eqtl_depth + 1
#         if sig_v4cordinate_dict[position].mpileup > 0:
#             number_eqtl_mpileup = number_eqtl_mpileup + 1
#         # print(sig_v4cordinate_dict[position])
#
#
#     print ("total_number_of_SNPs_loci\t" + str(len(snps_dict)))
#     missing_number = len(snps_dict) - len(sig_v4cordinate_dict)
#     print ("number_of_missing_SNPs_in_V4_genotype\t" + str(missing_number))
#     print ("number_of_non-missing_SNPs_in_V4_genotype\t" + str(len(sig_v4cordinate_dict)))
#     print ("totalDepth\t" + str(totalDepth))
#     print ("totalMpileup\t" + str(totalMpileup))
#     print ("totalChrLength\t" + str(len(seq)))
#     print ("number_snp_depth\t" + str(number_eqtl_depth))
#     print ("number_snp_mpileup\t" + str(number_eqtl_mpileup))
