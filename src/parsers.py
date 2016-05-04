import gzip
from numpy import *
import re
import itertools
import vcf

# convenience function to open coupled map/ped files
def plink_open(basename) :
    try :       mapfile = open(basename+'.map') 
    except :    mapfile = gzip.open(basename+'.map.gz')
    
    try :       pedfile = open(basename+'.ped')
    except :    pedfile = gzip.open(basename+'.ped.gz')

    return mapfile, pedfile

def map_parser(mapfile) :
    print '# Parsing MAP file :', mapfile.name
    mapraw = loadtxt (mapfile, dtype=int, usecols=(0,3))

    return sorted([(line[0], line[1]) for line in mapraw])

# parses a MAP file into a dictionary 
def map_dict(mapfile) : 
    print '# Parsing MAP file :', mapfile.name
    mapraw = loadtxt (mapfile, dtype=int, usecols=(0,3))

    mapdict = { i:[] for i in range(1,13)}
    for i in range(len(mapraw)) :
        mapdict[mapraw[i][0]].append([mapraw[i][1], i])

    # make sure it is sorted by position
    for i in range(1,13) :
        mapdict[i]=array(sorted(mapdict[i], key=lambda x:x[0]) ).T
        #mapdict[i]=array(mapdict[i] ).T
    
    return mapdict

# returns SNP indices and positions for a given (closed) interval
def map_find_loci(mapdict, sid, interval) :
    pos = mapdict[sid][0]
    indices = where((pos>=interval[0]) & (pos<=interval[1])) [0]
    lociSNPind = [mapdict[sid][1][i] for i in indices ]
    lociSNPpos = [mapdict[sid][0][i] for i in indices ]
    return lociSNPind, lociSNPpos


# reads ped file line by line and returns a string of SNPs 
# returns only homozygous SNPs (including '0')
# returns ' ' for heretozygous 
# optionally returns only SNPs from the given index

from itertools import izip,islice

# extract homozygous snps from a ped line
# this includes missing (0) calls
# replaces hetero snps with 

def pedline_homo(line) :
    # split line
    line = line.strip().split()
    
    # first column is sample name
    name = line[0]

    # split line into two strands
    # (be careful, data may not be phased)
    strand1 = islice(line, 6, None, 2)
    strand2 = islice(line, 7, None, 2)
    
    seq = zeros((len(line)-6) /2, dtype='|S1')
   
    for i, (alelle1, alelle2) in enumerate(izip(strand1,strand2)) :
        # homozgyous
        if alelle1 == alelle2 : seq[i] = alelle1 
        # heterozygous
        else : seq[i] = ' '

    return name, seq

# iterates over a ped file returning homozygous 

def ped_iterator(pedfile, index=None) :

    for line in pedfile :
        # read line
        # first column is sample name

        name, seq = pedline_homo(line)
        
        # optionally filter for index
        if index is not None :
            seq = seq[index]
        
        yield name, seq

# get the line for a given cultivar
# much faster than iterator because does not split every line
# identical output

def ped_find_cultivar( pedfile , cultivar ) :
    print '# Parsing PED file :', pedfile.name
    
    for line in pedfile :
        if line[:len(cultivar)] == cultivar :
            
            name, seq = pedline_homo(line)

            pedfile.seek(0)
            return name, seq

# parses a ped file 
# puts it in a '|S1' numpy array
# returns a name list and the array
# this is as memory efficient as it gets 
# still requires the whole thing fits into memory though

def ped_parser_homo(pedfile, nrows=3023, index=None) :
    print '# Parsing PED file :', pedfile.name
    
    snps=None
    names = []
    row=0

    for name, seq in ped_iterator(pedfile, index=index) :
        
        # initialize matrix after reading first line
        # this way you know nsnps
        if  snps is None:
            nsnps = (len(seq))
            snps = zeros([nrows,nsnps], dtype='S1')
        
        names.append(name)
        
        snps[row] = seq 
        row +=1
        if row==nrows : break
    
    print '# {} SNPs in {} cultivars.'.format(nsnps, nrows )
    return names,snps


# better vcf parser using pyvcf

def vcf_homo_snp (vcffile, minQual=50) : 
    reader = vcf.Reader(vcffile) 

    nsamples = len(reader.samples)

    assert nsamples == 1, 'function only implemented for single sample vcf files so far\n{}'.format(vcffile.name) 

    sample = reader.samples[0]

    positions = []
    sequence = ''
    for record in reader : 
        if not record.is_monomorphic and record.is_snp and record.QUAL>minQual :
            g1,g2 = record.genotype(sample)['GT'].split('/')

            if g1==g2 :
                g1=int(g1)
                assert g1 == 1, '{} genotype at {},{}'.format(record.genotype(sample)['GT'], record.CHROM, record.POS)

                chrno =int(re.search("\d+",record.CHROM).group())
                positions.append((chrno, record.POS))
                
                # weird struct here but ok.
                # maybe polish later

                assert len(record.ALT[0].sequence) == 1, 'ALT is weird.'
                
                sequence += record.ALT[0].sequence

    return positions, sequence

