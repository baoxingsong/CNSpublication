"""
Created: 25/12/2017  00:40:39
Author: Baoxing Song
Email: songbaoxing168@163.com
This source code is partially adopted from https://github.com/bvilhjal/mixmogam/
And the R source code of emma is referred

modified: 31 Oct 2019
By Baoxing Song try to take the CNS present/absent as genotypic variant for GWAS analysis using a fixed model
about the f test implemented here https://towardsdatascience.com/fisher-test-for-regression-analysis-1e1687867259
"""

import numpy as np
from numpy import linalg
from scipy import stats
import warnings
import re

class LinearModel:
    """
    A class for linear models
    """
    def __init__(self, Y=None, dtype='float'):
        """
        The fixed effects should be a list of fixed effect lists (SNPs)
        """
        self.n = len(Y)
        self.Y = np.matrix(Y, dtype=dtype)
        self.Y.shape = (self.n, 1)
        self.X = np.matrix(np.ones((self.n, 1), dtype=dtype))  # The intercept
        self.p = 1

    def updatey(self, Y, dtype='float'):
        assert (self.n == len(Y))
        self.Y = np.matrix(Y, dtype=dtype)
        self.Y.shape = (self.n, 1)
        # self.X = np.matrix(np.ones((self.n, 1), dtype=dtype))  # The intercept
        # self.p = 1
        # self.beta_est = None

    def add_factor(self, x, lin_depend_thres=1e-4):
        """
        Adds an explanatory variable to the X matrix.
        """
        new_x = np.array(x)
        new_x.shape = len(x)
        if lin_depend_thres>0.0:
            # Checking whether this new cofactor in linearly independent.
            (beta, rss, rank, sigma) = linalg.lstsq(self.X, new_x)
            if float(rss) < lin_depend_thres:
                warnings.warn(
                    'A factor was found to be linearly dependent on the factors already in the X matrix.  Hence skipping it!')
                return False
        new_x.shape = (self.n, 1)
        self.X = np.hstack([self.X, new_x])
        self.p += 1
        return True

    def get_estimates_fix_model(self, xs=None):
        X = np.hstack([self.X, xs])
        q = X.shape[1]  # number of fixed effects
        n = self.n  # number of individuls
        p = n - q
        (beta_est, rss, rank, sigma) = linalg.lstsq(X, self.Y)
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X, self.Y)
        if len(rss) == 0:
            return float(2.0)
        if len(h0_rss) == 0:
            return float(2.0)
        if h0_rss == rss:
            return float(2.0)
        f_stat = (h0_rss / rss - 1) * p / xs.shape[1]
        p_val = stats.f.sf(f_stat, (xs.shape[1]), p)
        if len(p_val) == 0:
            return float(2.0)
        # print(xs)
        # print("h0_rss:" + str(h0_rss) + " rss:" + str(rss) + " p:" + str(p) + " xs.shape[1]:" + str(xs.shape[1]) + " f_stat:" + str(f_stat) + " p_val:" + str(p_val))
        p_val = float(p_val)
        return p_val

def parse_cns_genotype_file(filename, absent_threshold=0.2, present_threshold=0.8, maf_count=5, min_total_count = 17):  # 25 peers and 3 pcs
    """
        read the cns length and coverage file and return a inter genotype table
    """
    individs = []
    with open(filename) as f:
        line = f.readline()
        l = list(map(str.strip, line.split()))
        for i in range(2, len(l)):
            individs.append(l[i])
    num_individs = len(individs)
    all_snps = {'chrs': [], 'positions': [], 'snps': []}
    with open(filename) as f:
        for line_i, line in enumerate(f):
            if line_i > 0:
                l = list(map(str.strip, line.split()))
                m = re.search('^(\d+):(\d+)\-(\d+)', l[0])
                if (m != None):
                    chrom = int(m.group(1))
                    snp = np.zeros(num_individs, dtype='float')
                    cns_length = float(l[1])
                    an = 0 # absent count
                    pn = 0 # present count
                    for i in range(2, num_individs+2, 1):
                        nt = float(l[i])/float(cns_length)
                        if nt <= absent_threshold:
                            snp[i-2] = 0.0
                            an = an + 1
                        elif nt >=present_threshold:
                            snp[i-2] = 1.0
                            pn = pn + 1
                        else:
                            snp[i-2] = 3.0 #missing
                    if( an >= maf_count and pn >= maf_count and (an+pn)>=min_total_count ):
                        all_snps['chrs'].append(chrom)
                        all_snps['positions'].append(int(m.group(2)))
                        all_snps['snps'].append(snp)
    return all_snps, individs

#read genotype data
genopype, individs = parse_cns_genotype_file("../../encodeCNSasGenotype/all.CNS.genotype.txt")
# get the genotypic variants taxa id to line number map. line number could be queried from the genotypic variants

import sys
maf_count = 5
min_total_count = 17
f = open(sys.argv[1], 'rU')
line_index = 0
phen_individs = {} # key is the taxa id value is the the index
phen_dict = {}  #key is the phenotype id
for line in f:
    if 0 == line_index:
        l = line.split()
        for i in range(0, len(l)):
            phen_individs[l[i]] = i
    else:
        d = []
        l = line.split()
        for i in range(1, len(l)):
            d.append(float(l[i]))
        pid = l[0]
        phen_dict[pid] = d
    line_index = line_index + 1

f.close()
#read phenotype data end

# get the PCS taxa id to line number map. line number could be queried from the genotypic variants
pcs = np.mat(np.loadtxt("../282set_PCs.txt", skiprows=1, usecols = (2,3,4))).astype("float")
pcs_id = {}
f1 = open("../282set_PCs.txt", 'rU')
i = 0
for line in f1:
    if i >0:
        l = list(map(str.strip, line.split()))
        pcs_id[l[0]] = i-1 # taxa id and the index of taxa in the matrix
    i = i+1

f1.close()
#pca reading end


lin_depend_thres=1e-4
for snp_index in range(len(genopype['snps'])):
    a_certain_gene = ""
    for gene in phen_dict:
        a_certain_gene = gene
        break
    p = phen_dict[a_certain_gene]
    this_pcs = []
    this_snp=[]
    this_phenotype = []
    used_index = []
    an = 0
    pn = 0
    for i in range(len(genopype['snps'][snp_index])):
        if (genopype['snps'][snp_index][i] != 3.0) and (individs[i] in pcs_id) and (individs[i] in phen_individs):
            used_index.append(i)
            if genopype['snps'][snp_index][i] ==0:
                an = an + 1
            else:
                pn = pn + 1
            this_snp.append(genopype['snps'][snp_index][i])
            this_pcs.append( pcs[pcs_id[individs[i]]] )
            this_phenotype.append(p[phen_individs[individs[i]]])

    if (an >= maf_count and pn >= maf_count and (an + pn) >= min_total_count):
        lmm = LinearModel(this_phenotype)
        gene_id = 0
        for gene in phen_dict:
            p = phen_dict[gene]
            this_phenotype = []
            for i in used_index:
                this_phenotype.append(p[phen_individs[individs[i]]])
            lmm.updatey(this_phenotype)
            if 0 == gene_id:
                this_pcs = np.mat(np.array(this_pcs))

                for i in range(np.size(this_pcs,1)):
                    lmm.add_factor(this_pcs[:, i], lin_depend_thres=lin_depend_thres)
            p_val = lmm.get_estimates_fix_model(xs=np.matrix(this_snp).T)
            gene_id = gene_id + 1
            if p_val <= 1.0:
                print( str(genopype['chrs'][snp_index]) + " " + str(genopype['positions'][snp_index]) + " " + str(p_val) + " " + gene)
