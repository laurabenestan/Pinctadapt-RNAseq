# RNA-seq using a reference genome

An integrated worklow based on genome-mapping and DE gene assessment to conduct RNA-seq data analyses in Colosse


**WARNING**

The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.


## Documentation


### Part 2: SNP analysis:

For this section, you will need to modify the cluster parameters and the $WORKDIR variable as well as the path to all sofwares according to your own specificities.


#### 1. Split NCigarStrings

Under construction


#### 2.Mark PCR Duplicates

Under construction


#### 3. Freebayes SNP calling:

We used Freebayes-parallel to perform SNP calling, with the following parameters: n-best-alleles=4, minCov=10 et minMapQ=30
For this, run:

```shell
01_freebayes-parallel.sh
```


#### 4. Filtering the VCF:

We used VCFlib, VCFtools and BCFtools to filter the VCf output of freebayes as follow:
* Keep only SNPs pattern (without complex events)
* A minimum depth (DP) of 10 reads per locus per genotype within an individual (under that, the genotype is transformed to "NA").
* A minor allele frequency (MAF) of at least 10% in the sampleset
* Less than 15% missing data per individual (over that, the individual is removed)
* Less than 10% missing data (miss) at a locus (over that, the locus is removed)
* Only loci that are biallelic
 

To filter the output (DP,MAF,miss), run one by one (waiting for the previous to be finished):

```shell
qsub 04_vcffilter.sh
qsub 05_vcftools.sh
qsub 06_BCFtools.sh
```

To remove individuals with too much missing data, we need to get the statistics from each individual.
For this run:

```shell
qsub 07a_vcftools_missing_Ind.sh
```
 
Then look at the output and list the individuals with more than 15% missing genotypes:
Modify the `07b_vcftools_missing_Ind.sh` file accordingly (list individuals to be removed), then run:

```shell
qsub 07b_vcftools_missing_Ind.sh
```
 
To remove complex SNPs (haplotype like genorype calling from Freebayes), run:

 
```shell
08_remove_complex_snps.sh
```


#### 5. Impute missing genotypes:

Now that the VCF is filtered, we have individuals with less than 15% missing data, and loci with less than 10% missing data.
We will use beagle5 to perform genotype imputation for the dataset.
For this, run:

```shell
qsub 09_beagle5_imputation.sh
```

 
#### 6. Perform eQTL analysis:

We will use XXXXXXX package to perform eQTL analysis.
To obtain the genotype table and the SNPlocation file necessary to run the analysis, we will use the VCF file produced by Beagle5. Run the script:
 
```shell
qsub 10_beagleVCF2eQTL.sh
```

Analysis under contruction........

 
 
#### 7. Perform outlier analysis:

To obtain candidate SNPs that differentiate the two populations, we used Bayescan to screen the loci.
We used default parameters with the exception of the prior odds (PO) of the neutral model, that we set to 500 (10 default). 
This parameter indicates our skepticism about the possibility that a given locus is under selection. 

For example, a PO = 10 indicates that we think the neutral model is 10 times more likely than the model with selection. The larger the PO, the more conservative the test of selection is. In principle, this parameter will not influence the results much if we have a large amount of data including many populations and many individuals per population. However, very frequently the number of populations or sample sizes are limited (e.g. less than 20 populations) so we do need to pay attention to this parameter.

The PO value that should be used depends on how many loci are included in the data set. If there are less than 1000 loci, then PO = 10 is reasonable, but with more loci (say between 1000 and 10000) PO = 100 or larger is a better choice. With millions of markers, as is the case in GWAS, values as large as 10000 may be necessary.

Here we have a bit more than 27000 loci, so we fixed PO at 300.

Bayescan requires its own input format. To build the correct input file, we need allelic frequencies within each population, the ploidy, and the number of alleles.

First, we will split the beagle5 imputed VCF file into two VCF files, one for each population:
For this, modify the following file according to your populations, and run:

```shell
qsub 11_singlePOP_VCF.sh
```

Then, run one after the other:
 
```shell
qsub 14_make_bayescan_input.sh
qsub 15_bayescan.sh
```

The script `14_make_bayescan_input.sh`will call the `make_bayescan_input.R`R code. This rscript requires two packages to be installed: *vcfR* and *plyr*.



#### 8. Visualize results of outlier analysis:

The R script `Bayescan_results_analysis.R`will use the ".sel" and "Fst" outputs from freebayes to visualize the significant outliers detected by bayescan. For this, we need the "coda" package to be installed.

```shell
qsub 16_get_bayescan_results.sh
```

This script will generate two plots for outlier visualisation, make a list of bayescan locus names with corresponding VCF locus names, and a table listing the significant outlier SNPs with their statistics.



#### 9. Prepare files for SNPeff analysis on the significant SNPs

For this step, we need to create a new VCf file that contains only the significant SNPs.

```shell
qsub 17_Filter_main_VCF_on_singificant_SNPs.sh
```

Now that we have the new VCF,we can run SNPeff on it:

```shell
18_SNPeff_on_significant.sh
```



#### EXTRA: Get statistics on the populations (heterozygosity and nucleotide diversity):

For this, we use VCFtools and the two populations specific VCF files we created above.
Run:
 
```shell
qsub 12_VCFtools_stats_popVCF.sh
```
This produces pi diversity and heterozygosity files for each populations.



