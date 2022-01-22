#!/usr/bin/env python3
'''
Classify transcripts in possible SQANTI3 SC if deleted from GTF
Scipt for modifying GTF
Author: Jorge Martínez Tomás
'''

import os
import argparse
import time

def readTarget(f_name):
    # Leer el fichero con los transcripts ID a eliminar
    trans = []
    refTrans = []
    refGene = []
    f_in = open(f_name, 'r')
    
    for line in f_in:
        line = line.split()
        trans.append(line[0])
        refTrans.append(line[2])
        refGene.append(line[3])
    f_in.close()

    return trans, refTrans, refGene

def getGeneID(line):
    line_split = line.split()
    gene_id = line_split[line_split.index('gene_id') + 1]
    gene_id = gene_id.replace(';', '')
    gene_id = gene_id.replace('"', '')
    
    return gene_id

def getTransID(line):
    line_split = line.split()
    trans_id = line_split[line_split.index('transcript_id') + 1]
    trans_id = trans_id.replace(';', "")
    trans_id = trans_id.replace('"', "")

    return trans_id

def check4modif(info, trans, refTrans, refGene, untouchable, deleted):
    to_write = []
    
    for line in info:
        line_split = line.split()
        feature = line_split[2]
        
        if feature == 'gene':
            to_write.append(line)
        else:
            trans_id = getTransID(line)
            if trans_id in trans:
                if trans_id in untouchable:
                    to_write.append(line)
                else:
                    # Se elimina este transcrito y se guarda su referencia
                    deleted.add(trans_id)
                    idx = trans.index(trans_id)
                    untouchable.append(refTrans[idx])
            else:
                to_write.append(line)
    if len(to_write) > 1:
        return to_write
    else:
        return None
            

def modifyGTF(f_name_in, f_name_out, trans, refTrans, refGene):
    f_out = open(f_name_out, 'w')
    
    gene_info = []
    untouchable = []
    deleted = set()
    prev_gene = str()
    prev_trans = str()
    write_curr_trans = True

    with open(f_name_in, 'r') as gtf_in:
        for line in gtf_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                gene_id = getGeneID(line)
                if not prev_gene:
                    prev_gene = gene_id

                if prev_gene == gene_id:
                    gene_info.append(line)

                else:
                    to_gtf = check4modif(gene_info, trans, refTrans, refGene, untouchable, deleted)
                    if to_gtf:
                        for i in to_gtf:
                            f_out.write(str(i))
                    gene_info = [line]
                    prev_gene = gene_id
    
    to_gtf = check4modif(gene_info, trans, refTrans, refGene, untouchable, deleted)
    if to_gtf:
        for i in to_gtf:
            f_out.write(str(i))

    f_out.close()

    print(len(deleted))    


def main():
    os.chdir('/home/jorge/Desktop/simulacion/modifgtf')
    f_delete = '/home/jorge/Desktop/ISM.filter.genconde.class.txt'
    f_gtf = '/home/jorge/Desktop/simulacion/getSC/gencode.v38.annotation.gtf'
    f_out = 'gencode.v38.annotation.ISM.del.gtf'

    trans, refT, refG = readTarget(f_delete)
    modifyGTF(f_gtf, f_out, trans, refT, refG)


if __name__ == '__main__':
    t = time.time()
    main()
    print(time.time() - t)