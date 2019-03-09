# Input 

- The  first two input files are:

 `m6a_promoter_dir`  : The position of SNP on m6a related gene promoter.

 `free_promoter_dir` : The position of SNP on gene promoter.

They should be **.bed** format ,require a tab character delimitation for the detail columns:

1. **chrom** - The name of the chromosome .
2. **chromStart** - The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
3. **chromEnd** - The ending position of the feature in the chromosome. The *chromEnd* base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as *chromStart=0, chromEnd=100*, and span the bases numbered 0-99.
4. **Entrez ID** - The Ensembl Gene ID of the gene which SNP was located. 
5. **strand** - Defines the strand. Either `.` (no strand) or `+` or `-`.

#### An example 

```
chr1    628911  631010  ENSG00000230021    -
chr1    815371  817470  ENSG00000177757    +
chr1    823138  825237  ENSG00000228794    +
chr1    919593  921692  ENSG00000223764    -
chr1    959210  961309  ENSG00000188976    -
```

- The third input file is `reference_genome`,it can be download on Ensembl or NCBI.

# Output

- The output file is `m6a_vs_free`  ,and includes three columns.

1. **tissue** - The tissue information of the result.
2. **description**   -  The description about m6a.
3. **meanOddsratio** 

#### An example

```
muscle  m6a     0.199786
muscle  m6a     0.406499
muscle  m6a     0.824236
muscle  m6a     0.854379
muscle  m6a     0.514041
```

