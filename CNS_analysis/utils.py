import re
def readFastaFile(fastaFile): # do not upcase the DNA sequence please
    fastas = {}
    name = ""
    seq = []
    with open(fastaFile) as f:
        for line in f:
            m = re.search('^>(\S+)', line)
            if (m != None):
                if (len(name) > 0) & (len(seq) > 0):
                    s = ''.join(seq)
                    s = re.sub("\\s", "", s)
                    fastas[name] = s
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            fastas[name] = s
    return fastas

class EQTLMarker:
    def __init__(self, chr, v3cordinate, v4cordinate, miniPvalue, significantCounts, depth, mpileup):
        self.chr = chr
        self.v3cordinate = v3cordinate
        self.v4cordinate = v4cordinate
        self.miniPvalue = miniPvalue
        self.significantCounts = significantCounts
        self.depth = depth
        self.mpileup = mpileup
    def __str__(self):
        return (self.chr + "\t" + str(self.v3cordinate) + "\t" + str(self.v4cordinate) + "\t" + str(self.miniPvalue) + "\t" + str(self.significantCounts) + "\t" + str(self.depth) + "\t" + str(self.mpileup))


class SNP:
    def __init__(self, chr, v3cordinate, v4cordinate, depth, mpileup):
        self.chr = chr
        self.v3cordinate = v3cordinate
        self.v4cordinate = v4cordinate
        self.depth = depth
        self.mpileup = mpileup
    def __str__(self):
        return (self.chr + "\t" + str(self.v3cordinate) + "\t" + str(self.v4cordinate) + "\t" + str(self.depth) + "\t" + str(self.mpileup))



# for parameter parsing begin
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        return True

def str2bool2(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        return False
# for parameter parsing end
