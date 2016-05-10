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

### Algorithm

This program will : 

1. Parse the *homozygous* SNPs from both the vcf and the map/ped dataset
2. Find those SNP positions that exist in both the unknown individual and the dataset
3. Create SNP-sequences based on those intersecting positions
4. Calculate and report the normalized hamming distances between the unknown individual and each sample in the dataset

Normalization is done by the size of intersection, that is the number of (homozygous) SNPs that exist both in the unknown vcf and the dataset, which is also reported in the output.  

### Data

To download various SNP-Seek datasets, one can make use of the makefile in the data folder. 

```
$ cd data
$ make NB-core_v4  	# for NB-core dataset
$ make 3krg_filt_snp_v4 # for (larger) 3krg_filt dataset
```

For detailed description of these datasets please refer to the [SNP-Seek paper](http://nar.oxfordjournals.org/content/43/D1/D1023.full) and the respective [download page](http://oryzasnp-atcg-irri-org.s3-website-ap-southeast-1.amazonaws.com/).
These downloads are provided for reproducibility purposes and validated by SHA checksums. 
One can substitute any other PLINK formatted (map/ped) dataset or convert a vcf-based dataset to map/ped format via [PLINK](https://www.cog-genomics.org/plink2/).
Please refer to the PLINK documentation for details on this process. 

The SNP data for the unknown individual is expected to be in the VCF format. 
In our case, this data was produced by the workflow [provided here](https://github.com/huangc/WGvarSNP).

Please note that for both pieces of data, we are only making use of homozygous SNPs. 

### Usage

#### Running

The python script identify.py is the backbone of this tool, it requires 2 positional arguments. It will take in 

1. the **map/ped base file path** and look for a %.map/%.ped or %.map.gz/%.ped.gz file pair in the specified path.  
2. Either
  * the **vcf file path** for an unknown individual
  * the **name** of an individual already in the dataset
3. Optionally, a genomic range

##### Comparing an unknown individual against the dataset

The minimal parameters to run are the plink dataset base path and the vcf file path, such as :

```
$ src/identify.py data/NB-core_v4 -vcf unknown.vcf 
```

Alternatively, one can also specify a range. The format is chromosome, beginning, end. For example :
 
```
$ src/identify.py data/NB-core_v4 -vcf unknown.vcf -range 4 100000 500000
```

will only consider the 100K to 500K region on the 4th chromosome. 
Note that the range is 1 based, beginning inclusive, end exclusive. 
But this will output to stdout, you can also specify an output file.


```
$ src/identify.py data/NB-core_v4 -vcf unknown.vcf -range 4 100000 500000 -out distances.dat

```
##### Comparing a sample from dataset against all others

identify.py can also run with a plink dataset and a sample name, in which case it will compare this sample against all others in the given dataset. 

```
$ src/identify.py data/NB-core_v4 -name B007
```

will do just that. 
Note that this will not work unless a sample with a given name exists in the dataset. 
Additionally, one can specify a range and an output file. 

```
$ src/identify.py data/NB-core_v4 -name B007-range 4 100000 500000 -out distances.dat


#### Output

One can simply ignore the lines starting with # as they are status logs and convenience messages
The non-comment lines consist of a cultivar name and a normalized hamming distance between the unknown individual and that cultivar. 
After the distance list, the program will report the exact input paths and number of snps for reference. 
Note that the hamming distances are normalized by the number of SNPs in the intersection of unknown cultivar and the database. 
The last 5 lines of the output will be the best five matches and associated hamming distances, reported for convenience.

#### Sorting

A quick way of removing the comments from the output file :

```
awk '$1 != "#" ' distances.dat > dist-nocomment.dat
```

A quick way of sorting the resulting a file would be :


```
sort -k 2 -n dist-nocomment.dat > dist-sorted.dat
```

And of course, these could be chained such as :

```
awk '$1 != "#" ' dist.dat | sort -k 2 -n > dist-sorted.dat
```

which will remove comments and sort the list by distance at once. 






