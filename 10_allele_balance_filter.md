# 10. Allele-Balance Assessment and Filtering

Following the removal of alleles longer than one base, the dataset contained **839,809 SNPs across 258 individuals**.

Allele balance (AB) was evaluated to identify SNPs exhibiting unequal support from reads associated with reference and alternate alleles in heterozygous genotype calls. Substantial allele imbalance can be associated with sequencing errors, mapping artifacts, or other complications in variant calling.

Because *Pocillopora acuta* can exhibit clonality and potentially variable ploidy, AB filtering was evaluated separately from Hardy–Weinberg equilibrium filtering.

## 10.1. Extract allele-balance values

```bash
bcftools query \
  -f '%INFO/AB\n' \
  09_maf001_maxmeanDP500_LEN1.vcf.gz \
  > 10_allele_balance_values.txt
```

The distribution was examined in R:

```r
library(ggplot2)

ab <- read.table(
  "10_allele_balance_values.txt",
  header = FALSE,
  col.names = "AB",
  na.strings = "."
)

summary(ab$AB)

ggplot(ab, aes(x = AB)) +
  geom_histogram(
    binwidth = 0.02,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = c(0.25, 0.75),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "Allele-Balance Distribution",
    subtitle = "258 individuals | 839,809 SNPs",
    x = "Allele balance (AB)",
    y = "Number of SNPs"
  ) +
  theme_classic()
```

![Allele-balance distribution before filtering](figures/10_allele_balance_distribution.png)

## 10.2. Evaluate the allele-balance filter

The selected ezRAD filtering criterion was:

```text
(AB > 0.25 AND AB < 0.75) OR AB < 0.01
```

This retains SNPs with intermediate allele-balance values and includes a special exception for values below 0.01.

The exception was retained to reproduce the selected ezRAD filtering criterion, although its biological justification for *P. acuta* remains uncertain. In particular, AB values near zero should not automatically be interpreted as evidence of reliable allele balance.

| AB category | SNP count | Percentage | Filter outcome |
|---|---:|---:|---|
| AB < 0.01 | 49,923 | 5.94% | Retained |
| 0.01 ≤ AB ≤ 0.25 | 204,345 | 24.33% | Removed |
| 0.25 < AB < 0.75 | 578,479 | 68.88% | Retained |
| AB ≥ 0.75 | 7,062 | 0.84% | Removed |
| **Total** | **839,809** | **100%** | |

Overall, the AB filter retained **628,402 SNPs (74.83%)** and removed **211,407 SNPs (25.17%)**.

## 10.3. Apply allele-balance filtering

```bash
bcftools view \
  -i '(INFO/AB>0.25 && INFO/AB<0.75) || INFO/AB<0.01' \
  -Oz \
  -o 10_allele_balance_filtered.vcf.gz \
  09_maf001_maxmeanDP500_LEN1.vcf.gz

bcftools index -t 10_allele_balance_filtered.vcf.gz
```

The retained SNP count was confirmed as **628,402**, with **258 individuals** remaining.

## 10.4. Filtering decision

Apply the selected ezRAD allele-balance criterion, retaining 628,402 SNPs for subsequent mapping-quality filtering.

**Final output:** `10_allele_balance_filtered.vcf.gz`

**Final SNP count:** **628,402**


#### CHANGE TO 0.2 and 0.8 as Maddie recommended

## Look at SNPs between 0.01 and 0.25 and then look at all the SNPs at the second peak and compare the individuals present 

**Final individual count:** **258**

**Next step:** Assess and filter mapping quality.
