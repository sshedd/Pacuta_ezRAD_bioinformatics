### 4: Assessing site missingness

The updated dataset contained **258 individuals and 8,933,146 SNPs**.

### Calculating Site Missingness

Site missingness was calculated using VCFtools:

```bash
vcftools \
  --gzvcf 07_initial_filtered_258ind.vcf.gz \
  --missing-site \
  --out 04_site_missingness_258ind
```

The resulting `04_site_missingness_258ind.lmiss` file contained the proportion of missing genotypes (`F_MISS`) at each SNP.

### Visualizing the Distribution

![Site missingness histogram after removing five individuals](figures/04_site_missingness_258ind.png)

| Maximum site missingness | Minimum completeness | SNPs retained | Percentage retained |
|---|---|---:|---:|
| 10% | 90% | 144,403 | 1.62% |
| 20% | 80% | 698,614 | 7.82% |
| **30%** | **70%** | **1,334,122** | **14.93%** |
| 40% | 60% | 1,894,520 | 21.21% |
| 50% | 50% | 2,446,666 | 27.39% |
| 60% | 40% | 2,977,251 | 33.33% |
| 70% | 30% | 3,518,930 | 39.39% |

### Filtering Decision

A minimum site-completeness threshold of **70%** (maximum site missingness of 30%) was selected for subsequent SNP filtering.

This threshold was chosen as a balance between retaining a substantial number of SNPs and reducing missing genotype data. It retains approximately **1.33 million SNPs (14.93%)** from the initial-filtered dataset. At 70%, the slope began to change, so it seemed like an appropriate threshold. 
