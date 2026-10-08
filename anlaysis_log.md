# *Pocillopora acuta* ezRAD Bioinformatics Analysis Log

## Project overview

This repository documents the bioinformatics processing of ezRAD sequencing data for *Pocillopora acuta*, from read processing and reference alignment through quality control, species identification, and preparation of SNP datasets for downstream population genomic analyses.

The workflow is being conducted on the Moneta Linux server.

**Primary analysis directory:**
`/home/sshedd/working/pacuta/`

**Conda environment:**
`ezRAD_env`

**Number of samples:** 273

**Nuclear reference genome:**
`Pocillopora_acuta_HIv1.assembly.fasta`

---

## 1. Initial dDocent analysis

The dDocent pipeline was used to process paired-end ezRAD sequencing reads and align them to the *Pocillopora acuta* nuclear reference genome.

**Working directory:**
`/home/sshedd/working/pacuta/dDocent_2`

The nuclear reference genome was copied into the dDocent directory and named `reference.fasta`.

Four blank samples and unassigned reads were excluded from the analysis.

### Input sequencing files

Original sequencing files follow the naming convention:

- `sample.R1.fq.gz`
- `sample.R2.fq.gz`

Processed sequencing files include:

- `sample.F.fq.gz`
- `sample.R.fq.gz`

### Troubleshooting notes

During the dDocent workflow, an error involving too many open files was encountered.

The file descriptor limit was increased using:

```bash
ulimit -n 1000000
```

A subsequent dDocent run entered the FreeBayes variant-calling stage, which was stopped so BAM processing and quality control could be completed first.

---

## 2. BAM filename corrections

Following a dDocent rerun, some BAM filenames contained a duplicated read-group suffix:

`sample-RG-RG.bam`

These files were renamed to:

`sample-RG.bam`

BAM indexes were regenerated, and outdated indexes were removed.

---

## 3. BAM mate-information correction

BAM mate information was corrected using samtools fixmate.

**Intermediate directory:**
`dDocent_fixmate_Pacuta`

A total of 273 corrected BAM files were generated.

Validation confirmed that properly paired reads had the expected paired-read structure.

The intermediate fixmate directory was deleted after successful validation of the downstream filtered BAM files.

---

## 4. Minimum template-length filtering

Alignments were filtered to retain records with an absolute template length greater than 99 bp.

The following filtering logic was applied to the corrected BAM files:

```bash
samtools view -h "$i" |
awk 'BEGIN {OFS="\t"} /^@/ {print; next} {if ($9 > 99 || $9 < -99) print}' |
samtools view -bS -o ../dDocent_min_tlen_fixmate_Pacuta/"$i"
```

Filtered BAM files were subsequently indexed with `samtools index`.

**Final filtered BAM directory:**

`/home/sshedd/working/pacuta/dDocent_min_tlen_fixmate_Pacuta`

**Final BAM count:** 273

### Validation example

For sample `103_PA01`:

- Total retained reads: 1,230,708
- Properly paired reads: 1,230,708
- Retained read pairs: 615,354
- Singletons: 0

These results describe the filtered BAM, not the original unfiltered alignment.

---

## 5. Updated mapped BED file

An updated `mapped.bed` file was generated from the filtered BAM files using bedtools.

The workflow involved:

1. Converting BAM alignments into BED format using `bedtools bamtobed`.
2. Combining the individual BED files.
3. Sorting genomic intervals.
4. Merging overlapping intervals.

**Final output:** `mapped.bed`

**Number of merged intervals:** 17,196

This BED file represents genomic regions covered by the retained alignments and will be used in subsequent analyses as appropriate.

---

## 6. Mapping quality assessment

Mapping quality was assessed using `samtools flagstat` and `samtools stats` across all 273 filtered BAM files.

### Software versions

- samtools 1.23.1
- htslib 1.23.1

### 6.1 Percentage mapped and properly paired

Statistics were generated using:

```bash
for i in *-RG.bam
do
    samtools flagstat "$i" > ./mapqc/samtools_flagstat/"$i".flagstat.txt
done
```

The average percentages were calculated across samples.

| Metric | Mean |
|---|---:|
| Reads mapped | 100.00% |
| Properly paired reads | 100.00% |

**Important:** These statistics were calculated from BAM files that had already undergone filtering. They should not be interpreted as the original mapping success rate.

### 6.2 Mapping quality scores (MAPQ)

Detailed mapping statistics were generated using:

```bash
for i in *-RG.bam
do
    samtools stats "$i" > ./mapqc/samtools_stats/"$i".stats.txt
done
```

All 273 statistics files were generated successfully.

The percentage of reads with MAPQ ≥ 30 was calculated for each sample using the MAPQ frequency tables.

```bash
echo -e "sample\t%_MAPQ30" > percent_MAPQ30.txt

for i in *.bam.stats.txt
do
    echo "$i"
    a=$(grep '^MAPQ' "$i" | head -n +29 | awk '{s+=$3}END{print s}')
    b=$(grep '^MAPQ' "$i" | tail -n -31 | awk '{s+=$3}END{print s}')
    c=$(awk -v a="$a" -v b="$b" 'BEGIN {printf "%.2f\n", 100 * b/(a+b)}')
    echo "$i" "$c" >> percent_MAPQ30.txt
done
```

The summary file contained 274 lines: one header and 273 sample results.

The mean percentage was calculated using:

```bash
awk 'NR>1 {sum += $2; n++} END {if (n>0) print sum/n}' percent_MAPQ30.txt
```

**Mean percentage of reads with MAPQ ≥ 30: 96.002%**

This exceeds the approximate 70% benchmark mentioned in the ezRAD workflow document.

### Mapping QC summary

| Metric | Result |
|---|---:|
| Number of samples | 273 |
| Mean percentage mapped | 100.00% |
| Mean percentage properly paired | 100.00% |
| Mean percentage MAPQ ≥ 30 | 96.002% |
| Merged mapped intervals | 17,196 |

---

## 7. Mitochondrial species identification

**Status: Planned**

Before nuclear SNP calling and filtering, samples will be examined at mitochondrial diagnostic SNP positions to help identify potential non-target species.

Planned approach:

1. Locate the mitochondrial reference genome.
2. Map processed ezRAD reads to the mitogenome.
3. Generate mitochondrial consensus sequences where sufficient coverage exists.
4. Examine diagnostic SNP positions in Geneious.
5. Record species assignments and unresolved samples.

Previous species-identification workflow:

[Species_ID_ezRAD](https://github.com/sshedd/Species_ID_ezRAD)

---

## 8. SNP calling and filtering

**Status: Not started**

SNP calling and filtering will proceed after the mitochondrial species-identification assessment.

---

## Reproducibility notes

- Long-running commands are executed within named `screen` sessions.
- Raw FASTQ and BAM files are retained on Moneta, not uploaded to GitHub.
- Analysis commands, parameter choices, troubleshooting steps, and numerical results are documented here.
- Future updates will be added chronologically as each stage is completed.
