import re,sys,os
import logging
from Bio.SeqIO import parse
from Bio.Seq import Seq,translate
from Bio.SeqIO import parse

seqs = {}
fid = open('yarrowia1.fa','w')

for ret in parse('/home/db/genomes/yarrowia_lipolytica/GCF_000002525.2_ASM252v1_genomic.fna','fasta'):
    seqs[ret.id] = str(ret.seq)
#print(len(seqs))

seq_s = set()  
gene_id = {} 

for rem in open('/home/db/genomes/yarrowia_lipolytica/cds.log').readlines():
    if rem[0] != '#':
        allinf = re.split('\t',rem)
#        print(allinf)

        if allinf[2] == 'CDS':
            if len(re.findall('GeneID:(.+?);',allinf[8])) != 0:
                geneid = re.findall('GeneID:(.+?);',allinf[8])[0]
		
                if geneid not in gene_id.keys():
                    gene_id[geneid] = 1
                    st1 = int(allinf[3]) - 1
                    st2 = int(allinf[4])
                    if allinf[6] == "-":
                        fid.write(geneid + '\t' + str(gene_id[geneid]) + '\t'+ allinf[3] +  '\t'+
									 allinf[4] + '\t'+ allinf[6] + '\t' + str(Seq(seqs[allinf[0]][st1:st2]).reverse_complement()) + '\n')
                    else:
                        fid.write(geneid + '\t' + str(gene_id[geneid]) + '\t'+ allinf[3] +  '\t' + allinf[4] + '\t'+ 
									allinf[6] + '\t' + seqs[allinf[0]][st1:st2] + '\n')
                else:
                    gene_id[geneid] += 1
                    st1 = int(allinf[3]) - 1
                    st2 = int(allinf[4])
                    if allinf[6] == "-":
                        fid.write(geneid + '\t' + str(gene_id[geneid]) + '\t'+ allinf[3] +  '\t'+ allinf[4] + '\t'+ 
								allinf[6] + '\t' + str(Seq(seqs[allinf[0]][st1:st2]).reverse_complement()) + '\n')
                    else:
                        fid.write(geneid + '\t' + str(gene_id[geneid]) + '\t'+ allinf[3] +  '\t' + allinf[4] + '\t'+ 
								allinf[6] + '\t' + seqs[allinf[0]][st1:st2] + '\n')

fid.close()

                    
