## Post-Individual-Removal MAC Reassessment

Following removal of five individuals with ≥90% missing genotypes and application of the 70% site-completeness threshold, the resulting VCF contained **258 individuals and 1,334,122 SNPs**.

Although a minor allele count (MAC) threshold of ≥2 had been applied during initial SNP filtering, removing individuals can change allele counts and cause previously polymorphic sites to become invariant or fall below the original MAC threshold.

### Checking for invariant sites

The number of SNPs with at least one observed minor allele was calculated using:

```bash
bcftools view \
  -c 1:minor \
  -H 09_site_completeness_70.vcf.gz | wc -l
```

**Result:** 1,334,010 SNPs retained at MAC ≥1.

Therefore, **112 sites** no longer contained an observed minor allele after individual removal.

A second diagnostic evaluated the original MAC ≥2 criterion:

```bash
bcftools view \
  -c 2:minor \
  -H 09_site_completeness_70.vcf.gz | wc -l
```

**Result:** 1,332,396 SNPs retained.

| Category | Number of SNPs |
|---|---:|
| SNPs before MAC reassessment | 1,334,122 |
| Invariant sites (MAC = 0) | 112 |
| Singleton minor alleles (MAC = 1) | 1,614 |
| Total sites failing MAC ≥2 | 1,726 |
| **SNPs satisfying MAC ≥2** | **1,332,396** |

Only 0.13% of SNPs failed the reassessed MAC ≥2 threshold.

### Reapplying the MAC filter

The MAC ≥2 criterion was reapplied to ensure that retained SNPs continued to satisfy the original filtering requirement.

```bash
bcftools view \
  -c 2:minor \
  -Oz \
  -o 10_site70_mac2.vcf.gz \
  09_site_completeness_70.vcf.gz

bcftools index -t 10_site70_mac2.vcf.gz
```

**Output:** `10_site70_mac2.vcf.gz`

**SNPs retained:** 1,332,396  
**Individuals retained:** 258

This VCF was used for subsequent SNP quality assessments.
