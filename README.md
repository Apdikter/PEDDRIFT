# *PEDDRIFT*
PEDDRIFT is an R tool for testing whether allele frequency differences between selection lines are consistent with genetic drift alone, or indicate selection pressure at specific loci.

## Overview

The method uses pedigree-based simulation to test the null hypothesis that observed genetic differentiation arose purely through drift, with no selection on the locus of interest. By simulating Mendelian inheritance through the actual pedigree structure, PEDDRIFT accounts for:

- Population structure and relatedness
- Sampling effects from limited genotyping
- Random genetic drift across generations
- Line-specific breeding schemes

## How it works

1. Takes observed genotype data from two selection lines
2. Simulates genetic inheritance through the pedigree under neutral evolution
3. Compares observed allele frequency differences to the null distribution
4. Identifies loci where differentiation exceeds what drift alone would produce

This pedigree-aware approach is more powerful than standard FST methods when detailed genealogical records are available, as it explicitly models the actual inheritance process rather than assuming random mating within subpopulations.

## Installation
```r
# Install from GitHub
devtools::install_github("Apdikter/PEDDRIFT")
```

## Quick Start
```r
library(PEDDRIFT)

# Load your pedigree data with genotypes
ped <- read.csv("pedigree_data.csv")

# Run PEDDRIFT with 1000 simulations
results <- peddrift(ped, nrep = 1000)

# View results
results$pvalues
```

## Input Requirements

### Required columns:
- `tag`: Unique animal identifier
- `sire`: Sire identifier (NA for founders)
- `dam`: Dam identifier (NA for founders)
- `line`: Selection line (must have exactly 2 lines coded as 1 and 2)
- `locus`: Genotype data (see formats below)

### Genotype formats:
- **Numeric**: 0, 1, 2 (homozygous reference, heterozygous, homozygous alternate)
- **Two-letter**: AA, AB, BB (or any two-letter combination)

### Important notes:
- **Pedigree must be topologically sorted**: All parents must appear before their offspring in the dataset
- **No duplicate animal IDs**: Each animal should appear exactly once
- **Founders**: Animals with missing parents (NA for both sire and dam) are treated as base population founders, with alleles sampled from the founder frequency distribution


### Parameters:

- **`pedigree`**: Data frame containing pedigree and genotype information
- **`afreq`**: Founder allele frequencies
  - `NULL` (default): Estimated as the mean frequency across lines
  - Single numeric value (0-1): Reference allele frequency (for biallelic loci)
  - Data frame: Custom frequencies with columns `allele` and `freq`
- **`nrep`**: Number of simulation replicates (default: 1000)
- **`rseed`**: Random seed for reproducibility (default: 0 = no seed)
- **`prune`**: Whether to prune pedigree to ancestors of genotyped animals (default: TRUE)
- **`allall`**: Whether to exclude simulations where alleles go extinct (default: FALSE)

### Output:

Returns a list containing:
- **`observed`**: Observed test statistics (chi-square, frequency difference, t-statistic)
- **`pvalues`**: P-values for each test statistic
- **`simulated`**: Data frame of simulated test statistics from all replicates
- **`summary`**: Summary information (number of lines, alleles, replicates, founder frequencies)


## Acknowledgements
Principal Scientist **Ken Dodds**, Bioeconomy Science Institute, AgResearch Group  
Principal Scientist **John McEwan**, Bioeconomy Science Institute, AgResearch Group
Summer Intern **Zak Morrison**, Bioeconomy Science Instititute, AgResearch Group

## References

A description of the theory and possible alternative methods is given in:  
Dodds, K G and McEwan, J C  (1997) Calculating exact probabilities of allele frequency differences in divergent selection lines. Proc Assoc Advmt Anim Breed Genet 12, 556-560. https://doi.org/10.13140/2.1.4568.7043


