## find best matching individual based on SNP fingerprint

### Aim 

This is a simple workflow to compare an individual described by SNPs in a given vcf file to a PLINK map/ped format SNP database to find the best matching individual. 
It was originally designed to match unknown samples to the reference cultivars from the [SNP-Seek database](http://nar.oxfordjournals.org/content/43/D1/D1023.full) based on [The 3000 Rice Genomes Project](http://gigascience.biomedcentral.com/articles/10.1186/2047-217X-3-7). 
The code is provided here as a general tool, it can be used to identify an indifidual given that one has:

  * SNP data for an unknown individual
  * SNP data for a reference population 

### Requirements

  * A working python (2) installation including libraries
    * numpy
    * pyvcf 
  * [SNP data based on 3krg project] (http://oryzasnp-atcg-irri-org.s3-website-ap-southeast-1.amazonaws.com/) is what we originally use this with and can be downloaded in map/ped format
    * The makefile in the data folder will donload and checksum these for you
  * [PLINK 1.9](https://www.cog-genomics.org/plink2/) can be used to convert any other vcf formatted SNP data to map/ped format required for the workflow

### Mechanism 

The python script identify.py will parse only the homozygous SNPs from the vcf and calculate the hamming distance for each cultivar in the database.
It will then report the normalized hamming distances for each sample in the database.
Normalization is done by the number of (homozygous) SNPs that exist both in the unknown vcf and the database, which is also reported in the output.  


### Data

To download the SNP-Seek database one can make use of the makefile in the data folder. 

```
$ cd data
$ make NB-core_v4  	# for NB-core dataset
$ make 3krg_filt_snp_v4 # for (larger) 3krg_filt dataset
```

For detailed description of these datasets please refer to the [SNP-Seek paper](http://nar.oxfordjournals.org/content/43/D1/D1023.full) and the respective [download page](http://oryzasnp-atcg-irri-org.s3-website-ap-southeast-1.amazonaws.com/)
These downloads are provided for reproducibility purposes and validated by SHA checksums. 
One can substitute any other PLINK formatted (map/ped) dataset or convert a vcf-based dataset to map/ped format via [PLINK](https://www.cog-genomics.org/plink2/).
Please refer to the PLINK documentation for details on this process. 

The SNP data for the unknown individual is expected to be in the VCF format. 
In our case, this data was produced by the workflow [provided here](https://github.com/huangc/WGvarSNP).

Please note that for both pieces of data, we are only making use of homozygous SNPs. 

### Usage

#### Running

The python script identify.py is the backbone of this tool, it has 2 required positional arguments. The identifier will take in 

1. the **map/ped base file path** and look for a %.map/%.ped or %.map.gz/%.ped.gz file pair in the specified path.  
2. the **unknown sample vcf file** path


```
$ src/identify.py data/NB-core_v4 unknown.vcf 
```

but this will output to stdout, you might want to pipe the output to a file


```
$ src/identify.py data/NB-core_v4 unknown.vcf > log.dat

```

#### Output

One can simply ignore the lines starting with # as they are status logs and convenience messages
The non-comment lines consist of a cultivar name and a normalized hamming distance between the unknown individual and that cultivar. 
After the distance list, the program will report the exact input paths and number of snps for reference. 
Note that the hamming distances are normalized by the number of SNPs in the intersection of unknown cultivar and the database. 
The last three lines of the output will be the best three matches and associated hamming distances, reported for convenience.

