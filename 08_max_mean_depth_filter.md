# 08. Maximum Mean Site Depth Filtering (500×)

Following global MAF filtering, a maximum mean site depth threshold of **500×** was applied to exclude exceptionally high-coverage SNPs that could potentially reflect repetitive genomic regions or mapping artifacts.

No additional minimum mean-depth threshold was applied because the mean site-depth distribution was relatively low and a minimum per-genotype depth of 3× had already been used during initial filtering.

## 08.1. Apply maximum mean site depth filter

```bash
vcftools \
  --gzvcf 07_global_maf001.vcf.gz \
  --max-meanDP 500 \
  --recode \
  --recode-INFO-all \
  --stdout \
  2> 08_max_meanDP500.log \
| bgzip -c > 08_maf001_maxmeanDP500.vcf.gz
```

## 08.2. Filtering results

| Metric | Result |
|---|---:|
| Input SNPs | 879,199 |
| SNPs removed | 13 |
| SNPs retained | **879,186** |
| Individuals retained | **258** |
| Maximum mean site depth | **500×** |

**Output:** `08_maf001_maxmeanDP500.vcf.gz`

The filter removed only **0.0015%** of SNPs, retaining nearly all markers while excluding the most extreme high-depth outliers.

## 08.3. Filtering decision

A maximum mean site depth threshold of 500× was retained because it removed exceptionally high-depth SNPs without substantially reducing the dataset.

**Final output:** `08_maf001_maxmeanDP500.vcf.gz`

**Next step:** Remove records containing reference or alternate alleles longer than one nucleotide.
