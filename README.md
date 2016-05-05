## find best matching individual based on SNP fingerprint

### Aim 

This is a simple workflow to compare a SNPs in a given vcf file to a PLINK format SNP database to find the best matching individual. 
It was originally designed to match unknown samples to the reference cultivars from the r3ksnpseek database. 
It can be generally used given that one has 
  * SNP data for a reference population 
  * SNP data for an unknown individual

### Requirements

  * a python (2) installation with libraries
    * numpy
    * pyvcf 
  * [SNP data based on 3krg project] (http://oryzasnp-atcg-irri-org.s3-website-ap-southeast-1.amazonaws.com/) is what we originally use this with and can be downloaded in map/ped format
    * The makefile in the data folder will donload and checksum these for you
  * [PLINK 1.9](https://www.cog-genomics.org/plink2/) can be used to convert any other vcf formatted SNP data to map/ped format required for the workflow

### Mechanism 

The python script identify.py will parse only the homozygous SNPs from the vcf and calculate the hamming distance for each cultivar in the database.
It will then report the normalized hamming distances for each sample in the database.
Normalization is done by the number of (homozygous) SNPs that exist both in the unknown vcf and the database, which is also reported in the output.  
The last three lines of the output will be the best three matches and associated hamming distances. 

### Usage

To download the SNPSeek database

```
$ cd data
$ make NB-core_v4  	# for core dataset
$ make 3krg_filt_snp_v4 # for (larger) filt dataset
```

To identify unknown sample described by unknown.vcf

```
$ src/identify.py data/NB-core_v4 unknown.vcf 
```
but this will print to stdout, you might want to pipe the output to a file
```
$ src/identify.py data/NB-core_v4 unknown.vcf > log.dat

```

