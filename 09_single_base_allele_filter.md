# 09. Allele-Length Filtering: Removal of Multinucleotide Alleles

Although the initial filtering pipeline selected SNPs using `bcftools -v snps`, additional inspection identified records with reference or alternate alleles longer than one nucleotide.

The VCF header contained an `INFO/LEN` annotation. However, allele lengths were assessed directly from the REF and ALT fields to ensure that only single-nucleotide alleles were retained.

## 09.1. Count records with multibase alleles

```bash
bcftools view \
  -i 'strlen(REF)!=1 || strlen(ALT)!=1' \
  -H 08_maf001_maxmeanDP500.vcf.gz \
| wc -l
```

This command counted records containing reference or alternate alleles longer than one nucleotide.

## 09.2. Remove multibase alleles

```bash
bcftools view \
  -e 'strlen(REF)!=1 || strlen(ALT)!=1' \
  -Oz \
  -o 09_maf001_maxmeanDP500_LEN1.vcf.gz \
  08_maf001_maxmeanDP500.vcf.gz

bcftools index -t 09_maf001_maxmeanDP500_LEN1.vcf.gz
```

## 09.3. Filtering results

| Metric | Result |
|---|---:|
| Input sites | 879,186 |
| Sites excluded | **39,377 (4.48%)** |
| SNPs retained | **839,809** |
| Individuals retained | **258** |
| Remaining records with allele length ≠ 1 | **0** |

**Output:** `09_maf001_maxmeanDP500_LEN1.vcf.gz`

## 09.4. Interpretation and decision

This step ensured that all retained variants were represented as single-nucleotide substitutions, improving compatibility with downstream population genomic and seascape genomic analyses.

**Filtering decision:** Exclude all records for which the REF or ALT allele length was not exactly one nucleotide.

**Final output:** `09_maf001_maxmeanDP500_LEN1.vcf.gz`

**Next step:** Assess and filter allele balance.
