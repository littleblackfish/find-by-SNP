#!/usr/bin/env python2
# coding: utf-8

from numpy import *
from parsers import *
from itertools import izip

#return hamming distance between strings
def hamdist(str1, str2):       
    assert len(str1) == len(str2) 
    return sum( ch1 != ch2 for ch1, ch2 in izip(str1, str2))


def pos_intersect(poslist1, poslist2, rng=None) :
   
    if poslist1 is poslist2 :
        intersection = set(poslist1)
    else :
        intersection = set(poslist1).intersection(set(poslist2))

    if rng is not None : 
        fullrange = set( [(rng[0], i) for i in range(rng[1], rng[2])])
        intersection = intersection.intersection(fullrange)


    assert len(intersection) != 0, 'Intersection is empty set.'
   
    index1=[]

    for i,pos in enumerate(poslist1) :
        if pos in intersection :
            index1.append(i)

    if poslist1 is poslist2 :
        return index1, intersection

    else :
        index2=[]
        for i,pos in enumerate(poslist2) :
            if pos in intersection : 
                index2.append(i)

        return index1, index2, intersection


if __name__ == '__main__' : 
    import argparse,sys,os

    parser = argparse.ArgumentParser(description='Program description')

    # plink file is mandatory
    parser.add_argument('plinkfname', help='basename for the map/ped file pair')

    # we can have an unknown vcf file
    parser.add_argument('-vcf', metavar='vcffile', help='vcf file for unknown individual')

    # alternatively, we can have a sample name
    parser.add_argument('-name', metavar='samplename', help='name of a sample of interest from map/ped file pair')

    # we can have a range
    parser.add_argument('-range', metavar=('CHR', 'BEG', 'END'),nargs=3, type=int, help='genomic range of interest')

    # we can optionally have an output file
    parser.add_argument('-out', metavar='outfile',type=argparse.FileType('w'), default=sys.stdout, help='optional output file (stdout by default)')

    args=parser.parse_args()

    if args.vcf is None and args.name is None :
          raise ValueError('There needs to be an vcf file or sample of interest in the plink file')

    # expand plink file base
    mapfile,pedfile = plink_open(args.plinkfname) 

    # this is the index for SNPs in the PED file
    plinkpos = map_parser(mapfile)

############### init complete, building reference ###################
    
    # if there is an unknown vcf
    if args.vcf is not None :         
        # strip vcf file to get homozygous snps only
        # positions are (chromosome, index) 
        with open(args.vcf) as vcffile :
            vcfpos, vcfseq = vcf_homo_snp(vcffile)
            unknownName = os.path.splitext(args.vcf)[0]


        # get intersection of SNPs in MAP and VCF files
        print '# Finding intersection...'
        vcfindex, mapindex, intersection = pos_intersect(vcfpos, plinkpos, args.range)


        #create reference sequence for the unknown cultivar
        refseq=''.join([vcfseq[i] for i in vcfindex])

            # if there is no unknown vcf, then we have a sample of interest in the dataset itself
    else :
        print '# Searching sample {} in plink dataset'.format(args.name)
        name, seq = ped_find_cultivar(pedfile, args.name)
        
        print '# Finding intersection...'
        mapindex, intersection = pos_intersect(plinkpos, plinkpos, args.range)


        # create reference sequence
        refseq=''.join( [seq[i] for i in mapindex] )

######### reference built, calculating distances ########

    # compare reference sequence with each cultivar in the PED file
    print "# Calculating hamming distances..."
    namelist=[]
    distlist=[]

    for cultName, cultSeq in ped_iterator(pedfile, index=mapindex) :
        namelist.append(cultName)
        dist = hamdist(refseq, cultSeq)/float(len(refseq))
        distlist.append(dist)
        args.out.write('{}\t{:.3f}\n'.format(cultName, dist))

    rank = argsort(distlist)

    # print some final stuff

    args.out.write('# PLINK dataset provded in {}\n'.format(args.plinkfname))
    if args.name is not None: args.out.write('# Sample of interest is {}\n'.format(args.name))
    if args.vcf is not None: 
        args.out.write('# Unknown cultivar provided in {}\n'.format(args.vcf))
        args.out.write('# {} SNPs in intersection of {} in unknown and {} in dataset.\n'.format(len(intersection),len(vcfpos),len(plinkpos)))
    
    # Print range
    if args.range is not None : 
        args.out.write('# Genomic range involved : chr {} [{},{})\n'.format(args.range[0], args.range[1], args.range[2]))

    # Print top 5
    args.out.write('# Closest 5 individuals\n')
    for i in range(5) :
        args.out.write('# {}\t{:.3f}\t{}\n'.format(i+1, distlist[rank[i]], namelist[rank[i]]))
    
    if args.out is not sys.stdout :
        args.out.close()
    print '# Done.'
