#!python
import re
import subprocess
import sys
from argparse import ArgumentParser
import sys
#read a fasta file and return a dictionary, the key is entry id and the value is the sequence in upcase
from utils import readFastaFile
from utils import EQTLMarker
import re
from utils import str2bool

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

    parser.add_argument("-e", "--eqtl",
                        dest="eqtl",
                        type=str,
                        default="",
                        help="eQtl result file")

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

    parser.add_argument("-p", "--pvalue",
                        dest="pvalue",
                        type=float,
                        default="0.00000001",
                        help="the pvalue threshold for eQTL significant")

    args = parser.parse_args()
    if args.genome == "":
        print("Error: please specify --genome", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.bam == "":
        print("Error: please specify --bam", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.eqtl == "":
        print("Error: please specify --eqtl", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.vcf == "":
        print("Error: please specify --vcf", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    reference_genome = readFastaFile(args.genome)
    print("reference genome reading done", file=sys.stderr)
    totalDepth = 0
    totalMpileup = 0

    reliable_genes = dict()
    with open("/media/bs674/2019junehackatho/eQTL/list") as f:
        for line in f:
            line = line.strip()
            reliable_genes[line] = 1

    sig_ids_dict = dict()
    chr = ""

    # read the eQTL result begin
    with open(args.eqtl) as f:
        for line in f:
            if not line.startswith("Trait"):
                elements = line.split('\t')
                if float(elements[6]) <= args.pvalue:
                    chr = elements[2]
                    if elements[0] in reliable_genes:
                        if elements[1] in sig_ids_dict:
                            sig_ids_dict[elements[1]].significantCounts = sig_ids_dict[elements[1]].significantCounts + 1
                            if sig_ids_dict[elements[1]].miniPvalue < float(elements[6]):
                                sig_ids_dict[elements[1]].miniPvalue = float(elements[6])
                        else:
                            e = EQTLMarker(chr, int(elements[3]), int(elements[3]), float(elements[6]), 1, 0, 0)
                            sig_ids_dict[elements[1]] = e
    # read the eQTL result end

    sig_v4cordinate_dict = dict()
    seq = reference_genome[chr]
    seq = re.sub("\\s", "", seq)
    seq = re.sub("-", "", seq)
    if args.mask:
        # read the VCF file begin
        with open(args.vcf) as f:
            for line in f:
                if line[0] is not '#':
                    elements = line.split('\t')
                    if elements[0] == chr:
                        elements2 = elements[2].split('-')
                        variant_id = "S" + elements2[0] + "_" + elements2[1]
                        if (variant_id in sig_ids_dict) and (reference_genome[chr][int(elements[1])-1] is not 'n'):
                            sig_v4cordinate_dict[elements[1]] = sig_ids_dict[variant_id]
                            sig_v4cordinate_dict[elements[1]].v4cordinate = int(elements[1])

        for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[2]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n'):
                    totalDepth = totalDepth + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].depth = int(elements2[2])

        for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[3]) > 0 and (reference_genome[chr][int(elements2[1])-1] is not 'n'):
                    totalMpileup = totalMpileup + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].mpileup = int(elements2[3])
        seq = re.sub("n", "", seq)
    else:
        # read the VCF file begin
        with open(args.vcf) as f:
            for line in f:
                if line[0] is not '#':
                    elements = line.split('\t')
                    if elements[0] == chr:
                        elements2 = elements[2].split('-')
                        variant_id = "S" + elements2[0] + "_" + elements2[1]
                        if (variant_id in sig_ids_dict) :
                            sig_v4cordinate_dict[elements[1]] = sig_ids_dict[variant_id]
                            sig_v4cordinate_dict[elements[1]].v4cordinate = int(elements[1])

        for line2 in subprocess.run(['samtools', 'depth', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[2]) > 0:
                    totalDepth = totalDepth + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].depth = int(elements2[2])

        for line2 in subprocess.run(['samtools', 'mpileup', '-r', chr, args.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout.split("\n"):
            if len(line2) > 0:
                elements2 = line2.split('\t')
                position = elements2[1]
                if int(elements2[3]) > 0:
                    totalMpileup = totalMpileup + 1
                if position in sig_v4cordinate_dict:
                    sig_v4cordinate_dict[position].mpileup = int(elements2[3])

    number_eqtl_depth = 0
    number_eqtl_mpileup = 0
    for position in sig_v4cordinate_dict:
        if sig_v4cordinate_dict[position].depth > 0:
            number_eqtl_depth = number_eqtl_depth + 1
        if sig_v4cordinate_dict[position].mpileup > 0:
            number_eqtl_mpileup = number_eqtl_mpileup + 1
        # print(sig_v4cordinate_dict[position])

    print ("total_number_of_eqtl_loci\t" + str(len(sig_ids_dict)))
    missing_number = len(sig_ids_dict) - len(sig_v4cordinate_dict)
    print ("number_of_missing_eqtl_loci_in_V4_genotype\t" + str(missing_number))
    print ("number_of_non-missing_eqtl_loci_in_V4_genotype\t" + str(len(sig_v4cordinate_dict)))
    print ("totalDepth\t" + str(totalDepth))
    print ("totalMpileup\t" + str(totalMpileup))
    print ("totalChrLength\t" + str(len(seq)))
    print ("number_eqtl_depth\t" + str(number_eqtl_depth))
    print ("number_eqtl_mpileup\t" + str(number_eqtl_mpileup))
