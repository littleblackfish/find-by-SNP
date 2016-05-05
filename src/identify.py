#!/usr/bin/env python2
# coding: utf-8

from numpy import *
from parsers import ped_iterator,map_parser,vcf_homo_snp,plink_open
from itertools import izip

#return hamming distance between strings
def hamdist(str1, str2):       
    assert len(str1) == len(str2) 
    return sum( ch1 != ch2 for ch1, ch2 in izip(str1, str2))


def pos_intersect(poslist1, poslist2) :
    
    intersection = sorted(set(poslist1).intersection(set(poslist2)))
    
    index1=[]
    index2=[]

    for i,pos in enumerate(poslist1) :
        if pos in intersection :
            index1.append(i)

    for i,pos in enumerate(poslist2) :
        if pos in intersection : 
            index2.append(i)

    return index1, index2, intersection

if __name__ == '__main__' : 

    import argparse,os

    parser = argparse.ArgumentParser() 
    parser.add_argument('plinkfname', help="plink file basename for the database")
    parser.add_argument('vcffname', help="vcf file for the unkown cultivar")
    parser.add_argument('-v', help="verbose mode, print as you calculate",
                    action="store_true")

    args=parser.parse_args()
            
    # strip vcf file to get homozygous snps only
    # positions are (chromosome, index) 
    with open(args.vcffname) as vcffile :
        vcfpos, vcfseq = vcf_homo_snp(vcffile)

    unknownName = os.path.splitext(args.vcffname)[0]

    # load plink MAP file 

    mapfile,pedfile = plink_open(args.plinkfname) 

    # this is the index for SNPs in the PED file
    plinkpos = map_parser(mapfile)

    # get intersection of SNPs in MAP and VCF files
    print '# Finding intersection...'
    vcfindex, mapindex, intersection = pos_intersect(vcfpos, plinkpos)

    #create reference sequence for the unknown cultivar
    vcfintseq=''
    for i in vcfindex :
        vcfintseq += vcfseq[i]

    # compare reference sequence with each cultivar in the PED file
    print "# Calculating hamming distances..."
    namelist=[]
    distlist=[]
    for cultName, cultSeq in ped_iterator(pedfile, index=mapindex) :
        namelist.append(cultName)
        dist = hamdist(vcfintseq, cultSeq)/float(len(intersection))
        distlist.append(dist)
        print '{}\t{:.3f}'.format(cultName, dist)

    rank = argsort(distlist)

    # write distance list
    with open(unknownName+'.dist', 'w') as f :
        print '# unknown cultivar provided in {}\n'.format(args.vcffname)
        print '# compared against database in {}\n'.format(args.plinkfname)
        print '# intersect\tunknown\tdatabase\n{}\t{}\t{}\n'.format(len(intersection),len(vcfpos),len(plinkpos))
        for i in range(3) :
            print '# {}\t{:.3f}\t{} \n'.format(i+1, distlist[rank[i]], namelist[rank[i]])
    
