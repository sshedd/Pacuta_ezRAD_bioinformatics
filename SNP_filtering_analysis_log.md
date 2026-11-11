SNP Filtering Analysis Log
================
2026-10-27

- [SNP Filtering Analysis Log](#snp-filtering-analysis-log)
  - [Project: *Pocillopora acuta* ezRAD Seascape
    Genomics](#project-pocillopora-acuta-ezrad-seascape-genomics)
  - [1. Objective](#1-objective)
  - [2. Starting Dataset](#2-starting-dataset)
  - [3. Software](#3-software)
  - [4. Working Directory](#4-working-directory)
  - [5. Filtering Strategy](#5-filtering-strategy)
  - [6. Initial SNP Filtering —
    BCFtools](#6-initial-snp-filtering--bcftools)
    - [Purpose](#purpose)
    - [Command](#command)
    - [Verification commands](#verification-commands)
  - [7. Next Filtering Stage —
    VCFtools](#7-next-filtering-stage--vcftools)
  - [8. Filtering Results Summary](#8-filtering-results-summary)
    - [Table 1. SNP filtering and quality-control
      summary](#table-1-snp-filtering-and-quality-control-summary)
    - [Table 2. Individual missingness
      distribution](#table-2-individual-missingness-distribution)
    - [Table 3. SNP completeness and potential retention
      thresholds](#table-3-snp-completeness-and-potential-retention-thresholds)
    - [Table 4. Mean SNP depth
      distribution](#table-4-mean-snp-depth-distribution)
    - [Table 5. Mean SNP depth among SNPs with ≥70%
      completeness](#table-5-mean-snp-depth-among-snps-with-70-completeness)
    - [Table 6. Output files generated](#table-6-output-files-generated)
  - [9. Additional Mean-Depth
    Assessment](#9-additional-mean-depth-assessment)
    - [9.1 Purpose](#91-purpose)
    - [Table 7. Upper-depth tail among SNPs with ≥70%
      completeness](#table-7-upper-depth-tail-among-snps-with-70-completeness)
  - [10. FreeBayes Mapping-Quality
    Assessment](#10-freebayes-mapping-quality-assessment)
    - [10.1 Purpose](#101-purpose)
    - [10.2 Mapping quality across all initially filtered
      SNPs](#102-mapping-quality-across-all-initially-filtered-snps)
    - [Table 8. Alternate-allele mapping quality across all
      SNPs](#table-8-alternate-allele-mapping-quality-across-all-snps)
    - [10.3 Mapping quality among SNPs with ≥70%
      completeness](#103-mapping-quality-among-snps-with-70-completeness)
    - [Table 9. Alternate-allele mapping quality among
      higher-completeness
      SNPs](#table-9-alternate-allele-mapping-quality-among-higher-completeness-snps)
    - [10.4 Comparing reference and alternate mapping
      quality](#104-comparing-reference-and-alternate-mapping-quality)
    - [Table 10. SNP retention under alternative mapping-quality
      thresholds](#table-10-snp-retention-under-alternative-mapping-quality-thresholds)
    - [10.5 Mapping-quality filtering
      decision](#105-mapping-quality-filtering-decision)
  - [11. Allele-Balance Assessment](#11-allele-balance-assessment)
    - [11.1 Purpose](#111-purpose)
    - [11.2 Allele-balance
      distribution](#112-allele-balance-distribution)
    - [Table 11. FreeBayes allele-balance
      distribution](#table-11-freebayes-allele-balance-distribution)
    - [11.3 Investigating SNPs with very low allele
      balance](#113-investigating-snps-with-very-low-allele-balance)
    - [Table 12. Heterozygous genotype counts among SNPs with AB
      \<0.1](#table-12-heterozygous-genotype-counts-among-snps-with-ab-01)
    - [11.4 Allele-balance filtering
      decision](#114-allele-balance-filtering-decision)
  - [12. Properly Paired Read
    Assessment](#12-properly-paired-read-assessment)
    - [12.1 Purpose](#121-purpose)
    - [12.2 Results](#122-results)
    - [Table 13. Properly paired read
      support](#table-13-properly-paired-read-support)
    - [12.3 Paired-read filtering
      decision](#123-paired-read-filtering-decision)
  - [13. Combined Candidate-Filter
    Assessment](#13-combined-candidate-filter-assessment)
    - [13.1 Purpose](#131-purpose)
    - [Table 14. Combined projected SNP
      retention](#table-14-combined-projected-snp-retention)
    - [13.2 Visualizing the combined filtering
      results](#132-visualizing-the-combined-filtering-results)
    - [13.3 Candidate filtering
      decision](#133-candidate-filtering-decision)
  - [14. Additional Visual Summaries](#14-additional-visual-summaries)
    - [Figure 6. Individual missingness across all 263
      samples](#figure-6-individual-missingness-across-all-263-samples)
    - [Figure 7. SNP retention by minimum sample
      completeness](#figure-7-snp-retention-by-minimum-sample-completeness)
  - [15. Updated SNP Filtering Decision
    Summary](#15-updated-snp-filtering-decision-summary)
  - [16. Reproducibility Notes](#16-reproducibility-notes)
  - [17. Next Planned Analyses](#17-next-planned-analyses)

# SNP Filtering Analysis Log

## Project: *Pocillopora acuta* ezRAD Seascape Genomics

**Started:** October 22, 2026  
**Status:** Initial SNP filtering planned / in progress  
**Repository:**
[Pacuta_ezRAD_bioinformatics](https://github.com/sshedd/Pacuta_ezRAD_bioinformatics)

------------------------------------------------------------------------

## 1. Objective

Filter raw SNPs generated by dDocent and FreeBayes to produce a
high-quality nuclear SNP dataset suitable for population genetic
analyses, species identification, technical-replicate assessment,
principal component analysis (PCA), and subsequent seascape genomics.

All 263 individuals will be retained during SNP filtering.

After SNP filtering, PCA and technical-replicate comparisons will be
used to evaluate species identification and sample quality. One suitable
individual per sampling site will subsequently be selected for
downstream seascape genomics analyses.

## 2. Starting Dataset

| Parameter              | Value                                    |
|------------------------|------------------------------------------|
| Organism               | *Pocillopora* spp., primarily *P. acuta* |
| Reference genome       | *Pocillopora acuta* HIv1                 |
| SNP caller             | FreeBayes through dDocent v2.9.8         |
| Number of individuals  | 263                                      |
| Raw variant sites      | 24,037,651                               |
| Input VCF              | `TotalRawSNPs.vcf`                       |
| Input VCF size         | Approximately 221 GB                     |
| SNP-calling completion | October 21, 2026                         |

The input for filtering is the original `TotalRawSNPs.vcf`, rather than
dDocent’s automatically filtered `Final.recode.vcf`.

## 3. Software

All analyses are performed in the `ezRAD_env` conda environment.

| Software             | Version                             |
|----------------------|-------------------------------------|
| VCFtools             | 0.1.17                              |
| BCFtools             | 1.23.1                              |
| vcflib (`vcffilter`) | Installed; version not yet recorded |

## 4. Working Directory

Original dDocent output:

``` text
/home/sshedd/working/pacuta/dDocent_min_tlen_fixmate_Pacuta/
```

SNP-filtering subdirectory:

``` text
/home/sshedd/working/pacuta/dDocent_min_tlen_fixmate_Pacuta/SNP_filtering/
```

The original raw VCF will be preserved. Filtering outputs will be stored
separately and named according to the filtering stage.

Available storage checked before filtering:

``` text
Filesystem: /dev/sda1
Total:      1.5 TB
Available:  1.3 TB
Usage:      7%
```

## 5. Filtering Strategy

The primary filtering workflow follows the supplied ezRAD analysis
guidelines, with additional consideration of filtering criteria used in
a published sea cucumber seascape genomics study.

Planned filtering stages:

1.  Initial technical filtering: biallelic SNPs, variant type, allele
    count, depth, and variant quality.
2.  Examine individual missingness and remove samples.
3.  Examine per-site missingness. 
4.  Apply additional SNP-level quality filters, including allele
    frequency and read-support criteria.
5.  Evaluate heterozygosity and potentially problematic loci.
6.  Reduce linked or closely spaced SNPs where appropriate.
7.  Produce a filtered dataset for PCA, species identification, and
    technical-replicate assessment.

**No individuals will be removed during the initial SNP-filtering
workflow.**

Because sampling sites contain only 1–3 individuals, population-level
MAF filters will not be applied using individual sampling sites as
populations.

Thresholds for site missingness, MAF, heterozygosity, and linkage will
be evaluated using the dataset and documented before application.

## 6. Initial SNP Filtering — BCFtools

### Purpose

Apply the initial SNP filters recommended in the ezRAD analysis document
while retaining all 263 individuals.

| Filter | Setting | Purpose |
|----|----|----|
| Maximum alleles | 2 | Retain biallelic variants |
| Variant type | SNPs | Exclude indels |
| Minimum non-reference allele count | 1 | Exclude invariant reference sites |
| Maximum mean genotype depth | 2,000 | Remove extremely high-depth sites |
| Require homozygous genotype | Not applied | Retain all-heterozygous sites for later evaluation |

### Command

``` bash
conda activate ezRAD_env

cd ~/working/pacuta/dDocent_min_tlen_fixmate_Pacuta

mkdir -p SNP_filtering
cd SNP_filtering

bcftools view \
    --max-alleles 2 \
    -v snps \
    -c 1:nref \
    -i 'MEAN(FORMAT/DP)<=2000' \
    ../TotalRawSNPs.vcf \
    -Oz \
    -o 01_biallelic_snps.vcf.gz
```

Compressed VCF output is used to reduce storage requirements.

### Verification commands

Run only after the filtering command completes successfully:

``` bash
bcftools index -t 01_biallelic_snps.vcf.gz

bcftools index -n 01_biallelic_snps.vcf.gz

bcftools query -l 01_biallelic_snps.vcf.gz | wc -l
```

## 7. Filtering Results Summary

### Table 1. SNP filtering and quality-control summary

| **Step** | **Analysis** | **Filters applied** | **SNPs retained** | **Individuals retained** | **Runtime** | **Status** |
|----|----|----|----|----|----|----|
| 0 | Raw FreeBayes variant calling | None | 24,037,651 | 263 | — | Complete |
| 1 | Initial BCFtools filtering | Biallelic SNPs; non-reference allele count ≥1; mean genotype depth ≤2,000 | 17,907,986 | 263 | Not recorded | Complete |
| 2 | Initial VCFtools filtering | `--minQ 20 --mac 2 --minDP 3` | 8,933,146 | 263 | 4,576 seconds (~76 min) | Complete |
| 3 | Individual missingness assessment | `--missing-indv` (diagnostic only) | 8,933,146 | 263 | 364 seconds (~6 min) | Complete |
| 4 | Per-site missingness assessment | `--missing-site` (diagnostic only) | 8,933,146 | 263 | 397 seconds (~7 min) | Complete |
| 5 | Per-site mean depth assessment | `--site-mean-depth` (diagnostic only) | 8,933,146 | 263 | 402 seconds (~7 min) | Complete |
| 5a | Mean depth among SNPs with ≥70% completeness | Combined site-missingness and mean-depth diagnostic files; no filtering applied | 1,256,391 SNPs assessed | 263 (unchanged) | Not recorded | Complete |
| 6 | Additional SNP filtering | MAF, site missingness, depth and other quality criteria to be determined | Pending | 263 planned | — | Not started |

### Step 1: Assess Individual Missingness

**Important:** All 263 coral samples are being retained throughout SNP
filtering. Individual selection will occur after PCA, species
identification checks, and technical-replicate assessment. Diagnostic
analyses do not modify the VCF or remove samples.

### Table 2. Individual missingness distribution

Calculated using VCFtools `--missing-indv` on 8,933,146 SNPs across 263
individuals.

| **Individual missingness** | **Number of individuals** | **Percentage** |
|----------------------------|--------------------------:|---------------:|
| Less than 10%              |                         0 |           0.0% |
| 10–30%                     |                         0 |           0.0% |
| 30–50%                     |                         5 |           1.9% |
| 50–70%                     |                       118 |          44.9% |
| 70–90%                     |                       135 |          51.3% |
| 90% or greater             |                         5 |           1.9% |
| **Total**                  |                   **263** |       **100%** |

## Individual Missingness Filtering: Removal of Five Individuals

### Rationale

Following a discussion with my advisor, we decided to remove individuals with greater than 90% missing genotype data before continuing site-level SNP filtering. Individuals with extremely high missingness contribute little usable genotype information and may influence the distribution of site missingness across the dataset.

Individual missingness was previously calculated using VCFtools from `02_initial_vcftools.vcf.gz`, which contained **263 individuals and 8,933,146 SNPs**. The resulting `03_individual_missingness.imiss` file was examined to identify individuals with at least 90% missing genotypes.

![Top 20 individuals ranked by missing genotype percentage](figures/08_top20_individual_missingness.png)

### Individuals identified for removal

Five individuals exceeded the missingness threshold:

| Individual | Proportion missing (`F_MISS`) | Percent missing |
|---|---:|---:|
| 189_PA01 | 0.914451 | 91.45% |
| 277_PA012 | 0.993308 | 99.33% |
| 332_PA012 | 0.912504 | 91.25% |
| 354_PA01 | 0.999562 | 99.96% |
| 372_PA01 | 0.950981 | 95.10% |

Two of these individuals (`277_PA012` and `332_PA012`) are technical replicate samples. Their exclusion will be accounted for in subsequent technical replicate analyses.

### Creating the exclusion list

The five sample IDs were saved in `07_high_missingness_samples.txt`:

```bash
printf '%s\n' \
  189_PA01 \
  277_PA012 \
  332_PA012 \
  354_PA01 \
  372_PA01 \
  > 07_high_missingness_samples.txt
```

### Removing individuals from the VCF

The five individuals were excluded directly from the initial-filtered VCF using `bcftools`. This avoids modifying the original FASTQ or BAM files or repeating the variant-calling pipeline.

```bash
bcftools view \
  -S ^07_high_missingness_samples.txt \
  -Oz \
  -o 07_initial_filtered_258ind.vcf.gz \
  02_initial_vcftools.vcf.gz \
  2> 07_remove_samples.log
```

The `-S ^` option instructs `bcftools` to exclude the individuals listed in the text file. The resulting VCF retains the remaining individuals and the original variant records.

The original 263-individual VCF was preserved for reproducibility.

### Input file integrity

An earlier attempt to remove the five individuals from `05_initial_vcftools_renamed.vcf.gz` failed with a BGZF read error, indicating a problem reading the compressed file.

The earlier file, `02_initial_vcftools.vcf.gz`, was subsequently tested:

```bash
bgzip -t 02_initial_vcftools.vcf.gz
echo $?
```

The test returned exit code `0`, indicating that the input passed the BGZF compression integrity check. Therefore, the sample exclusion was redirected to this earlier VCF.

### Verification of the filtered VCF

The expected output contains **258 individuals** and **8,933,146 SNP records**, because this step removes samples without applying additional site-level filters.

The following commands are used to verify the output:

```bash
bcftools index -t 07_initial_filtered_258ind.vcf.gz

bgzip -t 07_initial_filtered_258ind.vcf.gz
echo "Integrity test exit code: $?"

bcftools query -l 07_initial_filtered_258ind.vcf.gz | wc -l

bcftools index -n 07_initial_filtered_258ind.vcf.gz
```

These checks confirm file integrity, the number of retained individuals, and the number of variant records. **Final verification results should be recorded after the commands complete successfully.**

### 2: Assessing site missingness

## Site Missingness Analysis Following Individual Filtering

### Objective

Following the removal of five individuals with greater than 90% missing genotype data, site missingness was recalculated to evaluate the distribution of missing genotypes across SNPs and select an appropriate site-completeness threshold.

The updated dataset contained **258 individuals and 8,933,146 SNPs**.

### Calculating Site Missingness

Site missingness was calculated using VCFtools:

```bash
vcftools \
  --gzvcf 07_initial_filtered_258ind.vcf.gz \
  --missing-site \
  --out 08_site_missingness_258ind
```

The resulting `08_site_missingness_258ind.lmiss` file contained the proportion of missing genotypes (`F_MISS`) at each SNP.

VCFtools confirmed that 258 individuals and 8,933,146 sites were included in the analysis.

### Visualizing the Distribution

Site missingness was visualized using a histogram with a bin width of 0.02.

```r
library(ggplot2)

site_missingness <- read.table(
  "08_site_missingness_258ind.lmiss",
  header = TRUE
)

ggplot(site_missingness, aes(x = F_MISS)) +
  geom_histogram(
    binwidth = 0.02,
    boundary = 0,
    fill = "steelblue",
    color = "white"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Site Missingness After Removing 5 Individuals",
    subtitle = "258 individuals | 8,933,146 SNPs",
    x = "Proportion of individuals missing per SNP",
    y = "Number of SNPs"
  ) +
  theme_classic(base_size = 12)
```

![Site missingness histogram following removal of five individuals](figures/08_site_missingness_258ind.png)

**Figure:** Distribution of site missingness across 8,933,146 SNPs following removal of five individuals. A histogram bin width of 0.02 was selected to visualize the distribution.

The histogram showed a relatively flat distribution of SNP counts across intermediate missingness values, followed by a pronounced increase in SNP counts at high missingness values (approximately 0.85–1.00). No distinct break was evident that would independently justify a specific filtering threshold.

### Evaluating SNP Retention

The number of SNPs retained under several candidate maximum site-missingness thresholds was calculated in R:

```r
thresholds <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70)

retention <- data.frame(
  Max_missingness = thresholds,
  SNPs_retained = sapply(
    thresholds,
    function(x) sum(site_missingness$F_MISS <= x)
  )
)

retention$Percent_retained <- round(
  100 * retention$SNPs_retained / nrow(site_missingness),
  2
)

print(retention)
```

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

This threshold was chosen as a balance between retaining a substantial number of SNPs and reducing missing genotype data. It retains approximately **1.33 million SNPs (14.93%)** from the initial-filtered dataset.

Although the histogram did not reveal a clear natural breakpoint, 70% completeness provides a reasonable working threshold for downstream population genomic analyses.

**Decision:** Apply a minimum site completeness of 70% to the 258-individual dataset.

**Status:** Threshold selected; filtering will be performed in the next step.

### Table 5. Mean SNP depth among SNPs with ≥70% completeness

The site-missingness and mean-depth diagnostic files were compared to
assess sequencing depth among SNPs genotyped in at least 70% of the 263
individuals.

This analysis was diagnostic only; no additional SNPs or individuals
were removed.

| **Metric**                     |         **Result** |
|--------------------------------|-------------------:|
| SNPs meeting ≥70% completeness |          1,256,391 |
| Average mean SNP depth         |            12.227× |
| SNPs with mean depth below 3×  |                  0 |
| SNPs with mean depth 3–5×      |      3,148 (0.25%) |
| SNPs with mean depth ≥5×       | 1,253,243 (99.75%) |

**Interpretation:** SNPs meeting the provisional 70% completeness
threshold had substantially higher mean depth than the full dataset
(12.227× versus 4.753×). Nearly all SNPs in this subset (99.75%) had a
mean depth of at least 5×.

This indicates that low mean depth across the full dataset is strongly
associated with SNPs having higher missingness. Restricting analyses to
SNPs with greater completeness may improve the overall coverage profile
without requiring an additional minimum mean-depth filter.

**Filtering decision:** No additional minimum mean-depth threshold will
be applied at this stage. The provisional 70% completeness threshold
will be evaluated alongside minor allele frequency and other SNP-quality
criteria before final filtering.

### Table 6. Output files generated

| **Analysis** | **Output file** | **Purpose** |
|----|----|----|
| Initial BCFtools filtering | `01_biallelic_snps.vcf.gz` | Biallelic SNP dataset after initial filtering |
| Initial VCFtools filtering | `02_initial_vcftools.vcf.gz` | Dataset after quality, minor allele count and genotype depth filters |
| Individual missingness | `03_individual_missingness.imiss` | Missing genotype proportion for each sample |
| Site missingness | `04_site_missingness.lmiss` | Missing genotype proportion for each SNP |
| Mean SNP depth | `04_site_mean_depth.ldepth.mean` | Mean sequencing depth for each SNP |

## 9. Additional Mean-Depth Assessment

### 9.1 Purpose

The initial mean-depth assessment showed that SNPs with at least 70%
sample completeness generally had higher mean sequencing depth than the
full SNP dataset.

We therefore examined the upper tail of the mean-depth distribution to
determine whether unusually high-depth SNPs should be excluded.

Very high sequencing depth can sometimes indicate repetitive genomic
regions, duplicated sequences, or reads from multiple genomic locations
mapping to the same reference position. However, high depth alone does
not establish that a SNP is erroneous.

### Table 7. Upper-depth tail among SNPs with ≥70% completeness

| Mean SNP depth | Number of SNPs above threshold |
|----------------|-------------------------------:|
| \>50×          |                         11,109 |
| \>75×          |                          4,091 |
| \>100×         |                          1,884 |
| \>150×         |                            541 |
| \>200×         |                            167 |

**Interpretation:** Only 1,884 of the 1,256,391 SNPs meeting the
provisional 70% completeness threshold had mean depth greater than 100×,
representing approximately 0.15% of these SNPs.

**Provisional filtering decision:** Consider applying a maximum
mean-depth threshold of 100×.

**Reasoning:**

- Unusually high-depth SNPs may represent duplicated or repetitive
  genomic regions.
- The proposed threshold would remove only a small proportion of SNPs
  meeting the completeness requirement.
- A maximum-depth filter may reduce potentially problematic loci without
  substantially reducing the available SNP dataset.
- The 100× threshold is a pragmatic quality-control choice rather than a
  demonstrated biological boundary.

No additional minimum mean-depth threshold is proposed because the ≥70%
completeness subset already has substantially improved mean depth.

**Status:** Diagnostic complete. Maximum mean-depth filter not yet
applied.

## 10. FreeBayes Mapping-Quality Assessment

### 10.1 Purpose

Mapping quality indicates how confidently sequencing reads align to the
reference genome.

Two FreeBayes INFO annotations were examined:

- `MQM`: Mean mapping quality of reads supporting the alternate allele.
- `MQMR`: Mean mapping quality of reads supporting the reference allele.

Examining both metrics is important because reads supporting the
reference and alternate alleles may differ in mapping confidence.

### 10.2 Mapping quality across all initially filtered SNPs

The distribution of `MQM` was examined across all 8,933,146 SNPs in
`02_initial_vcftools.vcf.gz`.

### Table 8. Alternate-allele mapping quality across all SNPs

| MQM range | Number of SNPs | Percentage |
|-----------|---------------:|-----------:|
| \<20      |         96,582 |       1.1% |
| 20–30     |        298,545 |       3.3% |
| 30–40     |        482,484 |       5.4% |
| 40–50     |        896,227 |      10.0% |
| ≥50       |      7,159,308 |      80.1% |
| **Total** |  **8,933,146** |   **100%** |

**Interpretation:** Most SNPs had relatively high alternate-allele
mapping quality. Approximately 80.1% had `MQM ≥50`.

### 10.3 Mapping quality among SNPs with ≥70% completeness

We repeated the assessment for SNPs genotyped in at least 70% of the 263
individuals.

### Table 9. Alternate-allele mapping quality among higher-completeness SNPs

| MQM range | Number of SNPs |
|-----------|---------------:|
| \<20      |         10,308 |
| 20–30     |         36,603 |
| 30–40     |         39,913 |
| 40–50     |         96,047 |
| ≥50       |      1,073,520 |
| **Total** |  **1,256,391** |

Approximately 85.4% of these SNPs had `MQM ≥50`.

#### Figure 1. Alternate-allele mapping quality among SNPs with ≥70% completeness

``` mermaid
xychart-beta
    title "Alternate-allele mapping quality (MQM)"
    x-axis ["<20", "20-30", "30-40", "40-50", "50+"]
    y-axis "SNPs (millions)" 0 --> 1.2
    bar [0.010308, 0.036603, 0.039913, 0.096047, 1.073520]
```

**Figure interpretation:** Most SNPs meeting the completeness threshold
have high alternate-allele mapping quality. The higher-completeness
subset also has a larger proportion of high-MQM SNPs than the full
dataset.

### 10.4 Comparing reference and alternate mapping quality

We examined both `MQM` and `MQMR` among the 1,256,391 SNPs meeting the
provisional completeness threshold.

### Table 10. SNP retention under alternative mapping-quality thresholds

| Required mapping quality for both alleles | SNPs failing | SNPs passing |
|-------------------------------------------|-------------:|-------------:|
| ≥30                                       |       50,686 |    1,205,705 |
| ≥40                                       |       99,534 |    1,156,857 |
| ≥50                                       |      204,940 |    1,051,451 |

Additional diagnostic result:

**66,579 SNPs** had alternate-allele mapping quality more than 10 points
below reference-allele mapping quality.

This represents approximately 5.3% of SNPs meeting the completeness
threshold.

Such differences can indicate allele-specific mapping uncertainty,
although they do not automatically establish that a variant is
incorrect.

### 10.5 Mapping-quality filtering decision

**Provisional decision:** Require both `MQM ≥40` and `MQMR ≥40`.

**Reasoning:**

1.  Both reference and alternate alleles should have reasonably
    confident mapping support.
2.  Requiring both values to be at least 40 would retain approximately
    92.1% of SNPs meeting the completeness threshold.
3.  A threshold of 50 would remove substantially more SNPs.
4.  The proposed threshold may reduce false-positive SNPs associated
    with ambiguous read alignments.
5.  It balances quality control with retaining enough markers for PCA
    and species identification.

**Status:** Diagnostic complete. Mapping-quality filter not yet applied.

## 11. Allele-Balance Assessment

### 11.1 Purpose

Allele balance describes the relative representation of reference and
alternate reads among heterozygous observations.

For a diploid heterozygous genotype, approximately balanced read support
is generally expected, although sequencing depth, sampling variation,
and alignment bias can cause deviations.

FreeBayes `AB` is a site-level annotation, rather than a direct measure
of allele balance for each individual genotype.

### 11.2 Allele-balance distribution

The `AB` distribution was examined across all 8,933,146 initially
filtered SNPs.

### Table 11. FreeBayes allele-balance distribution

| AB range  | Number of SNPs |
|-----------|---------------:|
| \<0.1     |      1,312,114 |
| 0.1–0.2   |        420,636 |
| 0.2–0.3   |        983,509 |
| 0.3–0.4   |      1,567,457 |
| 0.4–0.6   |      3,991,528 |
| 0.6–0.8   |        632,086 |
| \>0.8     |         25,816 |
| **Total** |  **8,933,146** |

Approximately 44.7% of SNPs had AB between 0.4 and 0.6.

However, 1,312,114 SNPs had AB below 0.1, prompting further
investigation before choosing an allele-balance threshold.

### 11.3 Investigating SNPs with very low allele balance

Five SNPs with `AB <0.1` were inspected directly.

All five had `AB = 0` and no heterozygous genotype calls in the current
VCF. Instead, they contained combinations of homozygous-reference
genotypes (`0/0`), occasional homozygous-alternate genotypes (`1/1`),
and missing genotypes (`./.`).

We then counted heterozygous genotype calls across all SNPs with
`AB <0.1`.

### Table 12. Heterozygous genotype counts among SNPs with AB \<0.1

| Heterozygous calls per SNP | Number of SNPs | Percentage |
|----------------------------|---------------:|-----------:|
| Zero                       |      1,302,260 |     99.25% |
| 1–2                        |          8,343 |      0.64% |
| More than 2                |          1,511 |      0.12% |
| **Total**                  |  **1,312,114** |   **100%** |


**Figure interpretation:** Approximately 99.25% of SNPs with very low AB
had no heterozygous genotype calls in the current filtered VCF.

This is important because a SNP without heterozygous calls cannot be
meaningfully evaluated for read balance within its currently observed
heterozygotes.

An additional consideration is that FreeBayes calculated the original
INFO annotations before VCFtools masked genotypes with depth below 3×.
Therefore, the current genotype calls and original site-level AB
annotations may not reflect exactly the same observations.

### 11.4 Allele-balance filtering decision

**Decision: Do not apply a blanket site-level AB filter at this stage.**

**Reasoning:**

1.  Nearly all SNPs with AB below 0.1 have no heterozygous calls in the
    current VCF.
2.  Low AB therefore does not necessarily indicate poorly balanced
    heterozygous genotypes.
3.  SNPs without heterozygotes may still contain useful genetic
    variation.
4.  Such SNPs may be informative for distinguishing *Pocillopora*
    species during PCA.
5.  Applying a strict AB threshold could remove potentially valid
    markers unnecessarily.
6.  If allele balance becomes a concern later, genotype-level
    allele-balance assessment would be more informative.

**Status:** Diagnostic complete. No blanket AB filter recommended.

## 12. Properly Paired Read Assessment

### 12.1 Purpose

FreeBayes provides two annotations describing properly paired read
support:

- `PAIRED`: Fraction of alternate-allele observations supported by
  properly paired reads.
- `PAIREDR`: Fraction of reference-allele observations supported by
  properly paired reads.

Low values can sometimes indicate problematic alignments, although
paired-read support also depends on library preparation and read
processing.

### 12.2 Results

Both metrics were examined among SNPs meeting the provisional 70%
sample-completeness threshold.

### Table 13. Properly paired read support

| Diagnostic                            | Number of SNPs |
|---------------------------------------|---------------:|
| SNPs with ≥70% completeness           |      1,256,391 |
| Missing either paired-read annotation |              0 |
| Either allele paired fraction \<0.5   |              0 |
| Either allele paired fraction \<0.8   |              0 |
| Either allele paired fraction \<0.9   |              0 |
| Both paired fractions ≥0.9            |      1,256,391 |

All candidate SNPs had both `PAIRED ≥0.9` and `PAIREDR ≥0.9`.

**Interpretation:** Properly paired read support was consistently high
among the candidate SNPs.

Although this is reassuring, properly paired reads are not necessarily
uniquely or correctly mapped reads. Mapping quality remains an important
complementary measure.

### 12.3 Paired-read filtering decision

**Decision: Do not add a paired-read filter at this stage.**

**Reasoning:**

- All SNPs meeting the completeness threshold already have strong
  paired-read support.
- A threshold of 0.9 would remove no SNPs from the candidate dataset.
- Adding a filter that removes no SNPs would not improve the dataset.

**Status:** Diagnostic complete. No additional paired-read filter
recommended.

## 13. Combined Candidate-Filter Assessment

### 13.1 Purpose

Before creating another filtered VCF, we calculated the combined effects
of the proposed SNP-completeness, maximum-depth, and mapping-quality
thresholds.

This is important because SNPs may fail more than one quality criterion.
Adding removal counts from separate diagnostics could therefore
overestimate the total number of SNPs excluded.

The proposed criteria are:

1.  SNP completeness ≥70%.
2.  Mean SNP depth ≤100×.
3.  Alternate-allele mapping quality (`MQM`) ≥40.
4.  Reference-allele mapping quality (`MQMR`) ≥40.

### Table 14. Combined projected SNP retention

| Filtering stage                | SNPs remaining | SNPs removed at this stage |
|--------------------------------|---------------:|---------------------------:|
| Current initially filtered VCF |      8,933,146 |                          — |
| SNP completeness ≥70%          |      1,256,391 |                  7,676,755 |
| Plus mean SNP depth ≤100×      |      1,254,507 |                      1,884 |
| Plus MQM ≥40 and MQMR ≥40      |  **1,155,557** |                 **98,950** |

The depth and mapping-quality criteria together would exclude 100,834
SNPs after the completeness threshold.

Approximately 92.0% of SNPs meeting the 70% completeness threshold would
survive the additional proposed filters.

### 13.2 Visualizing the combined filtering results

#### Figure 3. SNP filtering workflow

``` mermaid
flowchart TD
    A["Raw FreeBayes variants<br/>24,037,651"]
    B["Biallelic SNP filtering<br/>17,907,986"]
    C["Initial VCFtools filtering<br/>8,933,146"]
    D["Quality diagnostics<br/>No additional SNPs removed"]
    E["Proposed: ≥70% completeness<br/>1,256,391"]
    F["Proposed: mean depth ≤100×<br/>1,254,507"]
    G["Proposed: MQM and MQMR ≥40<br/>1,155,557"]
    H["Future: species PCA and<br/>technical replicate assessment"]

    A --> B --> C --> D
    D -.-> E
    E --> F --> G --> H

    style A fill:#dbeafe,color:#17212b
    style B fill:#dbeafe,color:#17212b
    style C fill:#bbf7d0,color:#17212b
    style D fill:#e5e7eb,color:#17212b
    style E fill:#fef3c7,color:#17212b
    style F fill:#fef3c7,color:#17212b
    style G fill:#fef3c7,color:#17212b
    style H fill:#e5e7eb,color:#17212b
```

**Figure interpretation:** The first two filtering stages are complete.
The later thresholds have been evaluated using diagnostics but have not
yet been applied to generate a new VCF.

#### Figure 4. Projected SNP retention

``` mermaid
xychart-beta
    title "Projected SNP retention through candidate filtering"
    x-axis ["Initial VCF", "70% complete", "Depth <=100", "Both MQ >=40"]
    y-axis "SNPs (millions)" 0 --> 9
    bar [8.933146, 1.256391, 1.254507, 1.155557]
```

**Figure interpretation:** The 70% completeness threshold has the
largest effect on SNP retention. The proposed maximum-depth and
mapping-quality thresholds remove comparatively few additional SNPs.

#### Figure 5. Effect of depth and mapping-quality filtering after the completeness threshold

``` mermaid
xychart-beta
    title "Additional filtering after 70% completeness"
    x-axis ["70% complete", "Depth <=100", "Both MQ >=40"]
    y-axis "SNPs (millions)" 0 --> 1.4
    bar [1.256391, 1.254507, 1.155557]
```

**Figure interpretation:** After restricting SNPs to those present in at
least 70% of samples, the maximum-depth threshold removes only 1,884
SNPs. The mapping-quality criteria then remove another 98,950 SNPs.

### 13.3 Candidate filtering decision

**Provisional decision:** Generate a candidate filtered VCF using the
following thresholds:

| Criterion                        | Proposed threshold   |
|----------------------------------|----------------------|
| SNP completeness                 | ≥70%                 |
| Mean SNP depth                   | ≤100×                |
| Alternate-allele mapping quality | MQM ≥40              |
| Reference-allele mapping quality | MQMR ≥40             |
| Additional minimum mean depth    | None                 |
| Site-level allele balance        | No blanket filter    |
| Properly paired read support     | No additional filter |

**Reasoning:**

- The completeness threshold improves the overall depth profile while
  retaining more than 1.25 million SNPs.
- The maximum-depth filter removes a small number of unusually
  high-depth loci.
- The mapping-quality filters require reasonably confident read mapping
  for both alleles.
- The allele-balance diagnostic does not justify a blanket AB filter.
- The paired-read metrics already show strong support among all
  candidate SNPs.
- Approximately 1.16 million SNPs are projected to remain, providing a
  substantial candidate dataset for further evaluation.

**Important:** The projected SNP count is based on diagnostic
calculations. The candidate filtered VCF has not yet been generated, and
its actual SNP count must be verified.

## 14. Additional Visual Summaries

The following figures summarize earlier results from Sections 1–8. They
are included here so that the original tables and descriptions can
remain unchanged.

### Figure 6. Individual missingness across all 263 samples

``` mermaid
xychart-beta
    title "Individual genotype missingness"
    x-axis ["<30%", "30-50%", "50-70%", "70-90%", "90%+"]
    y-axis "Number of samples" 0 --> 150
    bar [0, 5, 118, 135, 5]
```

**Interpretation:** Most samples have high missingness in the current
VCF. Nevertheless, all 263 samples are being retained until
species-identification PCA and technical-replicate comparisons are
completed.

### Figure 7. SNP retention by minimum sample completeness

``` mermaid
xychart-beta
    title "SNP retention by sample completeness threshold"
    x-axis ["50%", "70%", "80%", "90%", "95%"]
    y-axis "SNPs (millions)" 0 --> 2.5
    bar [2.388252, 1.256391, 0.613206, 0.119552, 0.041854]
```

**Interpretation:** Increasing the minimum sample-completeness
requirement substantially reduces the number of retained SNPs.

The provisional 70% threshold represents a compromise between retaining
genetic markers and reducing missing data.

A stricter completeness threshold may be evaluated later as a
sensitivity analysis for PCA.

## 15. Updated SNP Filtering Decision Summary

| Filtering criterion | Status | Reason |
|----|----|----|
| Biallelic SNPs | Applied | Exclude indels and multiallelic variants |
| Minimum variant quality ≥20 | Applied | Remove lower-quality variant calls |
| Minimum minor allele count ≥2 | Applied initially | Exclude variants with very low observed minor allele counts |
| Minimum genotype depth ≥3× | Applied | Mask low-depth genotypes |
| Individual missingness | Diagnostic only | Retain all individuals for species PCA and replicate comparisons |
| SNP completeness ≥70% | Provisional | Balance SNP retention with reduced missingness |
| Additional minimum mean depth | Not recommended | Higher-completeness SNPs already show improved depth |
| Maximum mean depth ≤100× | Provisional | Remove a small number of unusually high-depth loci |
| MQM ≥40 and MQMR ≥40 | Provisional | Require good mapping quality for both alleles |
| Site-level allele balance | No blanket filter | Low AB is strongly associated with zero heterozygous calls |
| PAIRED and PAIREDR | No additional filter | All candidate SNPs have paired-read fractions ≥0.9 |
| Global MAF | Pending | Compare suitable thresholds before final downstream analyses |
| Heterozygosity filtering | Pending | Evaluate problematic loci without prematurely excluding species-informative SNPs |
| SNP thinning / linkage | Pending | Reduce non-independence where appropriate for PCA |
| Sample removal | Not permitted at this stage | Complete species identification and technical replicate assessment first |

## 16. Reproducibility Notes

The following points are important for interpreting and reproducing
these analyses:

- All diagnostics were performed using `02_initial_vcftools.vcf.gz`,
  unless otherwise stated.
- This VCF contains 8,933,146 SNPs and 263 individuals.
- The original `TotalRawSNPs.vcf` remains preserved.
- The candidate completeness, maximum-depth, and mapping-quality filters
  have not yet been applied.
- The site-missingness file `04_site_missingness.lmiss` contains six
  columns; `F_MISS` is column 6.
- The site-mean-depth file `04_site_mean_depth.ldepth.mean` contains
  four columns; `MEAN_DEPTH` is column 3.
- Combined diagnostics used chromosome and position checks to verify
  that corresponding rows from different files referred to the same SNP.
- FreeBayes INFO annotations were generated during variant calling and
  may not fully reflect subsequent genotype masking by VCFtools.
- All 263 coral samples remain retained.
- No population-level MAF filters will be applied using individual
  sampling sites as populations because sampling sites contain only 1–3
  individuals.
- Global MAF thresholds of 0.01 and 0.05 remain under consideration.
- The projected SNP counts are diagnostic estimates, not verified counts
  from a newly generated VCF.
- All Mermaid figures contain their numerical data directly in the
  Markdown code and do not require separate image files.

**Current status:** Initial SNP filtering and additional quality
diagnostics are complete. The next candidate filtering step has been
evaluated but not yet executed. All 263 individuals remain retained.

``` mermaid
pie showData
    title Individual missingness across 263 samples
    "30-50% missing" : 5
    "50-70% missing" : 118
    "70-90% missing" : 135
    "90% or more missing" : 5
```

## SNP Filtering Summary and Decision Log

### Overview

The initial dDocent variant-calling workflow produced **24,037,651 variant records across 263 samples**. SNP filtering was designed to retain reliable, informative variants for species identification, technical replicate comparisons, and subsequent seascape genomic analyses.

Because the dataset contains only 1–3 samples per sampling site, filtering decisions were based primarily on variant quality, sequencing depth, genotype completeness, and mapping reliability rather than population-specific allele frequencies or Hardy–Weinberg equilibrium.

**All 263 samples were retained throughout the filtering steps completed to date.** Sample selection will occur after SNP filtering, species identification using PCA, and technical replicate evaluation.

### Summary of filtering steps and results

| Step | Filter or analysis | SNPs retained | Reason for filter or analysis | Result and decision |
|---|---|---:|---|---|
| 0 | Original dDocent VCF (`TotalRawSNPs.vcf`) | 24,037,651 | Establish the starting dataset containing all called variants | 263 samples; baseline for filtering |
| 1 | Biallelic SNPs; exclude non-reference-only sites; mean site DP ≤2,000 | 17,907,986 | Retain simple SNPs and remove extreme-depth sites potentially associated with repetitive or problematic regions | 74.5% of original variants retained |
| 2 | QUAL ≥20; minor allele count ≥2; genotype DP ≥3 | 8,933,146 | Remove low-quality and singleton variants while masking individual genotype calls supported by fewer than three reads | 49.9% of Step 1 SNPs retained; 263 samples maintained |
| 3 | Individual missingness diagnostic | 8,933,146 | Assess variation in sequencing completeness among individuals before deciding on sample selection | 253 of 263 samples had ≥50% missing genotypes; no samples removed |
| 4 | Site missingness diagnostic; candidate ≥70% genotype presence | 1,256,391 projected | Identify loci with adequate genotype representation across samples for PCA and replicate comparisons | 14.1% of Step 2 SNPs meet the proposed threshold; filter not yet applied |
| 5 | Mean site depth diagnostic; candidate mean DP ≤100 | 1,254,507 projected | Exclude unusually high-depth sites that may represent repetitive regions, paralogs, or mapping artifacts | Removes 1,884 sites from the ≥70% presence subset; filter not yet applied |
| 6 | Mapping-quality diagnostic; candidate MQM ≥40 and MQMR ≥40 | 1,155,557 projected | Require reliable mapping of reads supporting both alternate and reference alleles | Combined with ≥70% presence and mean DP ≤100, retains approximately 1.16 million SNPs; filters not yet applied |
| 7 | Allele-balance diagnostic | No additional filtering | Evaluate whether allele balance could help identify unreliable SNPs | Most very low AB values occur at sites without heterozygotes; no blanket AB filter recommended |
| 8 | Paired-read diagnostic | No additional filtering | Evaluate whether poorly paired reads require additional filtering | All sites in the ≥70% presence subset had PAIRED and PAIREDR ≥0.9; no additional filter needed |

*Important: Steps 4–8 are diagnostic evaluations and proposed filtering decisions, not additional filtering commands already executed. The projected counts should not be described as final VCF results.*

### Key filtering decisions

**Biallelic SNPs:** Retaining biallelic variants simplifies downstream PCA and genotype-based analyses. Extreme mean depth was initially limited to 2,000 to exclude exceptionally high-coverage sites without prematurely applying a restrictive depth threshold.

**Variant quality and genotype depth:** A minimum site QUAL of 20 was selected as an initial quality threshold. A minor allele count of two removes singleton variants, while minimum genotype depth of three reduces reliance on very weak genotype calls. Low-depth genotypes are masked rather than removing entire samples.

**Site completeness:** A candidate threshold of 70% genotype presence was selected as a balance between SNP retention and missing data. The diagnostic results showed that approximately 1.26 million SNPs meet this threshold, providing a substantial marker set for further evaluation.

**Mean depth:** Sites with unusually high mean depth may reflect repetitive or duplicated genomic regions. A candidate upper threshold of 100 was supported by the observed depth distribution, where only 1,884 sites meeting the completeness threshold exceeded this value.

**Mapping quality:** Requiring both MQM and MQMR to be at least 40 provides a way to evaluate mapping reliability for both alleles. The combined diagnostic retained 1,155,557 candidate SNPs after the proposed completeness, depth, and mapping-quality requirements.

**Allele balance:** A blanket site-level AB filter was not selected because very low AB values were overwhelmingly associated with sites lacking heterozygotes. Applying such a filter could unnecessarily eliminate legitimate variants.

**Paired-read support:** Additional paired-read filtering was not selected because every site in the ≥70% presence subset already met the evaluated paired-read threshold.

**Sample retention:** All 263 samples remain in the dataset. Samples will not be removed until species identification and technical replicate comparisons have been completed.

### Overall retention

| Comparison | SNPs | Percentage |
|---|---:|---:|
| Original variant records | 24,037,651 | 100% |
| After initial biallelic SNP selection | 17,907,986 | 74.5% of original records |
| After initial VCFtools filtering | 8,933,146 | 37.2% of original records |
| Candidate sites with ≥70% genotype presence | 1,256,391 | 14.1% of Step 2 SNPs |
| Candidate sites meeting completeness, depth, and mapping-quality requirements | 1,155,557 | 12.9% of Step 2 SNPs |

The combined diagnostic suggests that approximately **1.16 million SNPs** could be retained under the proposed completeness, depth, and mapping-quality thresholds. This is a projected count, not yet a finalized filtered VCF.

### Next steps

1. Apply the selected site completeness, mean-depth, and mapping-quality filters to generate the next VCF.
2. Verify the resulting SNP count and confirm that all 263 samples remain present.
3. Evaluate additional filtering requirements, including linkage disequilibrium and genotype-level quality, as appropriate for the intended analyses.
4. Perform PCA to investigate species identification and genetic structure.
5. Compare technical replicates to assess genotyping consistency.
6. Select the appropriate sample representation for each sampling site before downstream seascape genomic analyses.

**Current status:** Initial filtering through Step 2 is complete. Site completeness, depth, mapping quality, allele balance, and paired-read diagnostics have been evaluated. The next combined SNP filtering step has not yet been applied.

