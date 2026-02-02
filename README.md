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




## Description
A description of the theory and possible alternative methods is given in:  
Dodds, K G and McEwan, J C  (1997) Calculating exact probabilities of allele frequency differences in divergent selection lines. Proc Assoc Advmt Anim Breed Genet 12, 556-560. https://doi.org/10.13140/2.1.4568.7043


## Installation
In R, `devtools::install_github("Apdikter/PEDDRIFT")`


## User Instructions
* User must supply a dataset containing pedigree information sorted such that all parents precede their progeny. For species which need to be at least one year old to have progeny, sorting on year of birth is usually sufficient. If not, a topological sort may be required.
* No animals should have duplicate entries or tags, and each animal specified as a parent should have its own entry.
* Animals with missing parents are treated as belonging to the base population, i.e., 'founders.' In this case, both alleles are sampled from the base population frequency provided in argument  `afreq`.

## Acknowledgements
Principal Scientist **Ken Dodds**, Bioeconomy Science Institute, AgResearch Group  
Principal Scientist **John McEwan**, Bioeconomy Science Institute, AgResearch Group
Summer Intern **Zak Morrison**, Bioeconomy Science Instititute, AgResearch Group
