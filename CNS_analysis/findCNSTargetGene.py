#!python
import re
import subprocess
import sys
from argparse import ArgumentParser
import sys
#read a fasta file and return a dictionary, the key is entry id and the value is the sequence in upcase
from utils import readFastaFile

class EQTLMarker:
    def __init__(self, chr, v3cordinate, v4cordinate):
        self.chr = chr
        self.v3cordinate = v3cordinate
        self.v4cordinate = v4cordinate
        self.targetGenes = []

if __name__ == '__main__':
    parser = ArgumentParser(description='infer the CNS and target gene using eQTL result,'
                                        'please input the vcf file and the eqtl for single chromosome only')
    parser.add_argument("-g", "--genome",
                        dest="genome",
                        type=str,
                        default="",
                        help="the masked reference genome file")

    parser.add_argument("-s", "--sam",
                        dest="sam",
                        type=str,
                        default="",
                        help="the output of and-CNS pipeline in bam format")

    parser.add_argument("-e", "--eqtl",
                        dest="eqtl",
                        type=str,
                        default="",
                        help="eQtl result file")

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

    if args.sam == "":
        print("Error: please specify --sam", file=sys.stderr)
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
    reliable_genes = dict()
    with open("/media/bs674/2019junehackatho/eQTL/list") as f:     # this is a list of genes with good quality eQTL
        for line in f:
            line = line.strip()
            reliable_genes[line] = 1

    sig_ids_dict = dict()

    # read the eQTL result begin
    with open(args.eqtl) as f:
        for line in f:
            if "Trait" not in line:
                elements = line.split('\t')
                if (elements[0] in reliable_genes) and (float(elements[6]) <= args.pvalue):
                    chr = elements[2]
                    if elements[1] not in sig_ids_dict:
                        e = EQTLMarker(chr, int(elements[3]), int(elements[3]))
                        sig_ids_dict[elements[1]] = e
                    sig_ids_dict[elements[1]].targetGenes.append(elements[0])
                    # print(elements[0] + "\t" + elements[1])
    # read the eQTL result end
    print("eQTL result reading done", file=sys.stderr)
    print("chr" + chr, file=sys.stderr)
    sig_chr_v4cordinate_dict = dict()
    # read the VCF file begin
    with open(args.vcf) as f:
        for line in f:
            if line[0] is not '#':
                elements = line.split('\t')
                if elements[0] not in sig_chr_v4cordinate_dict:
                    sig_chr_v4cordinate_dict[elements[0]] = dict()
                elements2 = elements[2].split('-')
                variant_id = "S" + elements2[0] + "_" + elements2[1]
                if variant_id in sig_ids_dict:
                    sig_chr_v4cordinate_dict[elements[0]][elements[1]] = sig_ids_dict[variant_id]
                    sig_chr_v4cordinate_dict[elements[0]][elements[1]].v4cordinate = int(elements[1])
                    # print(variant_id + "\t" + elements[2], file=sys.stderr)
    print("vcf file reading done", file=sys.stderr)

    v3_gene_id_to_v4_gene_id = dict()
    with open("/media/bs674/2t/testPan_cns/realStory/checkTheCNSAndTargetGeneUsingeQTLresult/gene_model_xref_v3.txt") as f:
        for line in f:
            elements = line.split("\t")
            if len(elements) > 11:
                v3_gene_id_to_v4_gene_id[elements[0]] = elements[10]

    with open(args.sam) as f:
        for line in f:
            line = line.strip()
            elements = line.split("\t")
            if len(elements) == 11:
                elements.append("V3")
                elements.append("V4")
            if len(elements) == 13:
                chr = elements[2]
                start = int(elements[3])
                end = start
                ms = re.findall('(\d+)[MDN=X]', elements[5])
                for mm in ms:
                    end += int(mm)
                v3Genes = dict()
                v3Genes.clear()
                if len(elements[11])>3:
                    v3Genes_array = elements[11].split(':')[1].split(";")
                    for e in v3Genes_array:
                        v3Genes[e] = 1
                if chr in sig_chr_v4cordinate_dict:
                    for position in range(start, end): # here end is not included, and the end value is larger by 1 than it should be
                        if str(position) in sig_chr_v4cordinate_dict[chr]:
                            for tg in sig_chr_v4cordinate_dict[chr][str(position)].targetGenes:
                                v3Genes[tg] = 1
                print('\t'.join(elements[:11]), end='')
                if len(v3Genes)>0:
                    print("\tV3:", end='')
                    output = ""
                    for gene in v3Genes:
                        output = output + gene+";"
                    output = output[:-1]
                    print(output, end='')
                    print("\tV4:", end='')
                    output = ""
                    for gene in v3Genes:
                        if gene in v3_gene_id_to_v4_gene_id:
                            output = output + v3_gene_id_to_v4_gene_id[gene] + ";"
                        else:
                            output = output + "-;"
                    output = output[:-1]
                    print(output, end='')
                print()
            else:
                print(line)
