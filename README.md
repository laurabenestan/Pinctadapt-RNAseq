# RNA-seq using a reference genome

An integrated worklow based on genome-mapping and DE gene assessment to conduct RNA-seq data analyses in Colosse


**WARNING**

The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.


## Documentation


### Part 2: SNP analysis:

#### 1. Freebayes SNP calling:
under construction

#### 2. Filtering the VCF:

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
Modify the `07b_vcftools_missing_Ind.sh` file accordingly, then run:

```shell
qsub 07b_vcftools_missing_Ind.sh
```

To remove complex SNPs (haplotype like genorype calling from Freebayes), run:

```shell
08_remove_complex_snps.sh
```


#### 3. Impute missing genotypes:
Now that the VCF is filtered, we have individuals with less than 15% missing data, and loci with less than 10% missing data.
We will use beagle5 to perform genotype imputation for the dataset.
For this, run:

```shell
qsub 09_beagle5_imputation.sh
```

#### 4. Perform eQTL analysis:







