# *PEDDRIFT*
PEDDRIFT is a tool for testing allele frequency differences between selection lines and allows for esitmation of the effect of genetic drift. The technique is based on simulating genetic inheritance under the null hypothesis (that establishment of, and selection within the lines has not placed any selection pressure on the locus of interest), based on the pedigree and sampling structure in the data. 


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
