## 2. Initial VCFtools

Following successful completion of the initial BCFtools stage, apply the following filters from the ezRAD analysis document:

- Minimum variant quality: `--minQ 20`
- Minimum minor allele count: `--mac 2`
- Minimum genotype depth: `--minDP 3`

### Apply initial VCFtools filters

```bash
vcftools \
  --gzvcf 01_initial_bcftools.vcf.gz \
  --minQ 20 \
  --mac 2 \
  --minDP 3 \
  --recode \
  --recode-INFO-all \
  --stdout | bgzip -c > 02_initial_vcftools.vcf.gz

bcftools index -t 02_initial_vcftools.vcf.gz
```

### Filtering results

| Metric | Result |
|---|---:|
| Individuals retained | 263 |
| SNPs retained | 8,933,146 |
| Minimum variant quality | 20 |
| Minimum minor allele count | 2 |
| Minimum genotype depth | 3 |

The `--minDP 3` option masks individual genotype calls with fewer than three supporting reads rather than removing the entire SNP solely because one individual has insufficient depth.

The resulting VCF, `02_initial_vcftools.vcf.gz`, was used for subsequent individual missingness filtering.
