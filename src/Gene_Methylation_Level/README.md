# Input 

- The  first two input files are:

` peak_m6A_dir`  : The directory of peak methylation level in each tissuer.

 `peak_gene_dir` : The peak and related gene in each tissue.

They should be **.bed** format ,require a tab character delimitation for the detail columns:

1. **chrom** - The name of the chromosome .
2. **chromStart** - The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
3. **chromEnd** - The ending position of the feature in the chromosome. The *chromEnd* base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as *chromStart=0, chromEnd=100*, and span the bases numbered 0-99.
4. **Entrez ID** - The Ensembl Gene ID of the gene which SNP was located. 

#### An example 

```
chr10   100151489       100152302       ENSG00000107566.13      
chr10   100230242       100230369       ENSG00000227492.1       
chr10   100248246       100248328       ENSG00000230224.1       
chr10   100248246       100248328       ENSG00000095485.17      
chr10   100347169       100347296       ENSG00000099194.5       
```

- The third input file is `len_file`. Length of longest insform each gene in gencode v27 annotation file.
  1. **Entrez ID** - The Ensembl Gene ID .
  2. **Length** - The length of longest insform each gene.

#### An example

```
ENSG00000242268.2       710
ENSG00000188026.12      3932
ENSG00000270112.3       4685
ENSG00000280143.1       5326
```



# Output

- The output file is in directory  `Gene_methylation`  ,and includes two columns.

1. **Entrez ID** - The Ensembl Gene ID .
2. **methylation level**   -  It is the  methylation level in each gene.

#### An example

```
ENSG00000242268.2       0.337470614085
ENSG00000058866.14      0.523552930624
ENSG00000270112.3       1.10716848858
ENSG00000280143.1       0.151822215359
```

