#!/usr/bin/env python
import collections
from Bio import SeqIO
import sys
'''
This a little hack from NLR out summarizer script by Phillip Bayer: https://gist.github.com/philippbayer/0052f5ad56121cd2252a1c5b90154ed1
Here it will print NLR classed of each output sequence defined by their motif organization.
The output will look like
SeqName NLR_class
scaffold10886_nlr_1 CNL
scaffold2236_nlr_1 CNL
scaffold2236_nlr_2 CNL
'''
# this the motifs table from https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-75
motifs_d = {1:'nb_arc_cnl_or_tnl',2:'nb_arc_cnl',3:'nb_arc_cnl_or_tnl',4:'nb_arc_cnl_or_tnl',5:'nb_arc_cnl_or_tnl',6:'nb_arc_cnl',7:'linker_cnl_or_tnl',8:'linker_cnl_or_tnl',9:'lrr_cnl_or_tnl',10:'nb_arc_cnl_or_tnl',11:'lrr_cnl_or_tnl',12:'nb_arc_cnl_or_tnl',13:'tir_tnl',14:'monocot',15:'tir_tnl',16:'prenb_cnl',17:'prenb_cnl',18:'tir_tnl',19:'lrr_cnl_or_tnl',20:'monocot'}

# this is my idea of class assignment
class_dict = {frozenset(['lrr_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'CNL',frozenset(['lrr_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCNL', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CNL',frozenset(['nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'monocot', 'prenb_cnl', 'nb_arc_cnl']):'CNL', frozenset(['monocot', 'nb_arc_cnl_or_tnl']):'N', frozenset(['linker_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']):'TCNL', frozenset(['monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'NL', frozenset(['nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TN',frozenset(['nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl']):'N', frozenset(['tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'NL', frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCN', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl']):'TN', frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CNL', frozenset(['monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'TCNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['nb_arc_cnl_or_tnl', 'prenb_cnl']):'CN', frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCNL', frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl']):'TN', frozenset(['linker_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl']):'N', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl']):'NL', frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl']):'NL', frozenset(['monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'TCN', frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']):'TNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']):'CNL', frozenset(['nb_arc_cnl_or_tnl']):'N', frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'TCNL', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl']):'CNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']):'CN', frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl']):'TCNL', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl']):'N', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CN', frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl']):'CN', frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']):'CNL'}

# a list or a string of input files - NLR-Parser txt output
files = sys.argv[1].split('\n')


from collections import defaultdict
all_types = ['TN', 'TCNL', 'NL', 'CN', 'N', 'TCN',  'TNL','CNL']
header = ['SeqName','NLR_class'] 
print('\t'.join(header))

for f in files:
    count_dict = {}
    for a in ['complete', 'complete (pseudogene)', 'partial', 'partial (pseudogene)']:
        count_dict[a] = {'TN':0, 'NL':0, 'CN':0, 'N':0, 'TCN':0, 'TCNL':0, 'TNL':0, 'CNL':0}


    for line in open(f):
        ll = line.rstrip().split('\t')
        seqname = ll[1]
        motifs = ll[-1].split(',')
        this_domains = set()
        # iterate ovr
        for m in motifs:
            m  = int(m)
            this_domains.add(motifs_d[m])
        this_domains = frozenset(this_domains)
        this_class = class_dict[this_domains]
        print(seqname,"\t",this_class)
#        print(this_class)

