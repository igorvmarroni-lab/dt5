# Data Sources

All simulation inputs are from publicly available sources.
**No individual-level genomic or clinical data are included in this repository.**
IRB waiver applies (simulation study; no human subjects data analysed).

\---

## 1\. UK Biobank Allele Frequencies — DT3 and DT5

|Field|Value|
|-|-|
|Source|IEU Open GWAS|
|Dataset ID|ukb-b-6906|
|N variants|153,914 PASS variants|
|Ti/Tv ratio|2.502|
|Access|https://gwas.mrcieu.ac.uk|
|Citation|Elsworth et al. (2020). *Int J Epidemiol* 49:1558–1568|

**How used:** Allele frequencies (AF) sampled to simulate genotype matrix
`G \~ Binomial(2, af\_j)` for N\_SNPS=500 variants per individual.
In the code, a Beta(1.5, 1.5) proxy is used when the real AF file is not
present; replace with the real ukb-b-6906 AF vector for exact reproduction.

**To obtain:** Use the `ieugwasr` R package:

```r
# install.packages("ieugwasr")
library(ieugwasr)
variants <- associations(id = "ukb-b-6906", ...)
```

\---

## 2\. BD GWAS Effect Sizes — DT3 Calibration

|Field|Value|
|-|-|
|Source|Psychiatric Genomics Consortium (PGC)|
|Publication|Mullins et al. (2021). *Nat Genet* 53:817–829|
|Calibration value|mean\|β\| = 0.077 (DT3 Half-Normal sigma = 0.044)|
|Note|Effect sizes are NOT per-SNP betas; DT3 samples from a calibrated distribution|

**Citation:** Mullins N, et al. Genome-wide association study of more than
40,000 bipolar disorder cases provides new insights into the underlying biology.
*Nature Genetics* 53, 817–829 (2021). https://doi.org/10.1038/s41588-021-00857-4

\---

## 3\. BD GWAS Summary Statistics — DT5 Calibration

|Field|Value|
|-|-|
|Source|Psychiatric Genomics Consortium (PGC3)|
|File|daner\_bip\_pgc3\_nm\_noukbiobank.gz|
|N|40,463 BD cases / 313,436 controls|
|Valid SNPs|7,119,536 (INFO ≥ 0.8)|
|Empirical mean\|β\||0.0168 (range: −0.49 to +0.42)|
|Access|https://pgc.unc.edu/for-researchers/download-results/|

**How used:** DT5 samples effect sizes from `Half-Normal(sigma=SIGMA\_DT5)` where
`SIGMA\_DT5` is derived so that `E\[|beta|] = 0.0168`. The per-SNP betas from the
summary statistics file were NOT used directly (no rsID crossmatch).

**Citation:** Same as above (Mullins et al. 2021).

\---

## 4\. MEPS Longitudinal Data — DT4

|Field|Value|
|-|-|
|Source|Agency for Healthcare Research and Quality (AHRQ)|
|Survey|Medical Expenditure Panel Survey (MEPS)|
|Years|1996–2023|
|N (full cohort)|878,642 person-year observations|
|Analytic cohort|N=28,052 non-cases aged 15–30 (PHQ2<3 at baseline)|
|Eligible cohort (DT4)|N=3,398 (K6SUM ≥ 7)|
|Access|https://meps.ahrq.gov/mepsweb/data\_stats/download\_data\_files.jsp|
|License|Public domain (US Government work)|

**How used:** Real MEPS data serve as the substrate for Plasmode simulation (DT4).
The analytic cohort must be constructed from raw MEPS panel files. Key variables:

* `PHQ2` (0–6): Patient Health Questionnaire-2 (sum of two items)
* `K6SUM` (0–24): Kessler Psychological Distress Scale
* Eligibility: age 15–30, PHQ2<3 at baseline, K6SUM≥7, valid PHQ2 at ≥2 rounds

**To reproduce:** Download MEPS longitudinal data files and run the cohort
construction script (not included; see MEPS documentation for variable names
across panel years).

**Critical note (from paper):** K6SUM≥7 is NOT equivalent to PRS-BD ≥75th
percentile. K6SUM captures current psychological distress (may be environmental).
PRS-BD captures inherited genetic liability. This non-equivalence is acknowledged
as a limitation and runs in a conservative direction.

\---

## Synthetic Substitutes (for code validation)

When real data files are not present, the scripts generate synthetic substitutes:

* **DT3/DT5:** UKB AF proxy via `Beta(1.5, 1.5)` truncated to MAF \[0.01, 0.49]
* **DT4:** Synthetic MEPS cohort preserving N=3,398 and 14.9% transition rate

**Synthetic data produce results for code validation only.**
They do not reproduce the paper's reported NNT/RR estimates.
Use real data sources above for exact reproduction.

