# SNP Filtering Summary — *Pocillopora acuta* ezRAD

This table summarizes the sequential SNP quality-control and filtering decisions made during the *Pocillopora acuta* ezRAD bioinformatics pipeline.

| Step | Analysis | Filtering decision / criteria | SNPs removed | SNPs retained | Individuals retained |
|---|---|---|---:|---:|---:|
| 00 | Raw FreeBayes variant calling | No filtering applied | — | **24,037,651** | 263 |
| 01 | Initial BCFtools filtering | Biallelic SNPs; non-reference allele count ≥1; mean genotype depth ≤2,000 | 6,129,665 | **17,907,986** | 263 |
| 02 | Initial VCFtools filtering | `--minQ 20 --mac 2 --minDP 3` | 8,974,840 | **8,933,146** | 263 |
| 03 | Individual missingness assessment | Assess missing genotype proportions; no SNP filtering | 0 | **8,933,146** | 263 |
| 04a | Remove high-missingness individuals | Exclude 5 individuals with high missingness | 0 | **8,933,146** | 258 |
| 04b | Site completeness filtering | Retain sites genotyped in ≥70% of remaining individuals | 7,599,024 | **1,334,122** | 258 |
| 05 | Reapply minor allele count | Global MAC ≥2 after site completeness filtering | 1,726 | **1,332,396** | 258 |
| 06 | Mean site depth assessment | Examine depth distribution; no SNP filtering | 0 | **1,332,396** | 258 |
| 07 | Global minor allele frequency | Global MAF ≥0.01 | 453,197 | **879,199** | 258 |
| 08 | Maximum mean site depth | Mean site depth ≤500× | 13 | **879,186** | 258 |
| 09 | Single-base allele filtering | REF and ALT allele lengths must equal 1 | 39,377 | **839,809** | 258 |
| 10 | Allele-balance filtering | `(AB > 0.25 AND AB < 0.75) OR AB < 0.01` | 211,407 | **628,402** | 258 |
| 11 | Mapping-quality filtering | MQM >30; MQMR >30; MQM/MQMR between 0.75 and 1.25 | 32,596 | **595,806** | 258 |
| 12 | Heterozygosity and HWE assessment | No HWE filter; heterozygosity cutoff not yet selected | 0 | **595,806** | 258 |

## Overall filtering outcome

- **Initial SNPs:** 24,037,651
- **Final SNPs retained:** 595,806
- **Total SNPs removed:** 23,441,845
- **Percentage of SNPs retained:** 2.48%
- **Percentage of SNPs removed:** 97.52%
- **Initial individuals:** 263
- **Final individuals:** 258
- **Individuals removed:** 5

## Interpretation

The greatest reductions in SNP number occurred during initial quality filtering and site completeness filtering.
