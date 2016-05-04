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
  * [[PLINK 1.9]](https://www.cog-genomics.org/plink2/) can be used to convert any other vcf formatted SNP data to map/ped format required for the workflow
