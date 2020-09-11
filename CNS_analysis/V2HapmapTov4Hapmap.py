#!python
import re
from argparse import ArgumentParser
import sys
#read a fasta file and return a dictionary, the key is entry id and the value is the sequence in upcase
from pyliftover import LiftOver
# https://pypi.org/project/pyliftover/
# please run it with python3

if __name__ == '__main__':
    parser = ArgumentParser(description='uplift hapmap file from V2 coordinate to V4')

    parser.add_argument("-m", "--hapmap",
                        dest="hapmap",
                        type=str,
                        default="",
                        help="hapmap file with V2 coordinate")

    parser.add_argument("-c", "--chain",
                        dest="chain",
                        type=str,
                        default="",
                        help="the chain file downloaded from: ftp://ftp.ensemblgenomes.org/pub/plants/release-44/assembly_chain/zea_mays/AGPv2_to_B73_RefGen_v4.chain.gz ")
    args = parser.parse_args()


    if args.hapmap == "":
        print("Error: please specify --hapmap", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.chain == "":
        print("Error: please specify --chain", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    lo = LiftOver(args.chain)

    # read the SNP hapmap begin
    with open(args.hapmap) as f:
        for line in f:
            line = line.rstrip()
            if line[0] is 'S':
                elements = line.split('\t')
                chr = elements[2]
                position = int(elements[3]) - 1 # the coordinates in the liftover tools are 0 based
                lf = lo.convert_coordinate(chr, position)
                if None != lf and len(lf) == 1:
                    newChr = lf[0][0]
                    newPosition = lf[0][1] + 1 # change the 0 based coordinate to 1 based coordinate
                    print(elements[0] + "\t" + elements[1] + "\t" + newChr + "\t" + str(newPosition), end = '')
                    for i in range(4, len(elements)):
                        print("\t" + elements[i], end = '')
                    print()
            else:
                print(line)
