# Input

- The  three input files are:

 `eQTL SNP` : The position of  eQTL SNPs, and will enrich it.

 `m6a_bed_dir`  : The position of m6a region.

 `con_bed_dir` : The control file of the m6a region.

They should be **.bed** format ,require a tab character delimitation for the detail columns:

1. **chrom** - The name of the chromosome .
2. **chromStart** - The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
3. **chromEnd** - The ending position of the feature in the chromosome. The *chromEnd* base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as *chromStart=0, chromEnd=100*, and span the bases numbered 0-99.

####  An example 

```
chr1	14641	15015
chr1	52041	52415
chr1	87571	87945
chr1	108889	109263
chr1	160127	160501
```



# Output

- The output file is `Fisher_exact_test_result`  ,and includes three columns.

1. **tissue** - The tissue information of the result.
2. **trait**   -  The trait of the tissue.
3. **meanOddsratio** 

#### An example

```
tissue	trait	meanOddsratio
brain	brain_amygdala	1.680870
brain	brain_anterior_cingulate_cortex_ba24	1.701867
brain	brain_caudate_basal_ganglia	1.660124
brain	brain_cerebellar_hemisphere	1.702149
```





