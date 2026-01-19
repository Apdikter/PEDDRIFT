# PEDDRIFT
PEDDRIFT tests for allele frequency differences between selection lines and allows for the effect of genetic drift. The technique is based on simulating genetic inheritance under the null hypothesis (that establishment of, and selection within the lines has not placed any selection pressure on the locus of interest), based on the pedigree and sampling structure in the data. 

A description of the theory and possible alternative methods is given in Dodds, K G and McEwan, J C  (1997) Calculating exact probabilities of allele frequency differences in divergent selection lines. Proc Assoc Advmt Anim Breed Genet 12, 556-560. https://doi.org/10.13140/2.1.4568.7043


* The user must supply a dataset containing pedigree information. This must be sorted so that parents precede progeny. For species which need to be at least one year old to have progeny,  sorting on year of birth or including year of birth as a prefix in the individual identification field, and using this as the sort criterion, is usually sufficient. Otherwise, a topological sort may be required.






* No animals should have duplicate entries or tags, and each animal specified as a parent should have its own entry. The program checks for duplicate entries only for those animals which are parents. In this case, or in the case where a specified parent does not have its own entry, a message will be issued and the program will stop after checking the remainder of the pedigree file.

* The program treats animals with missing parents as those that belong to the base population. If an animal has one parent unknown, that parent is treated as a member of the base population. In the latter case the program randomly samples one allele from the base population allele frequencies and one from the known parent, in the former case both alleles are sampled from the base population.



CONSIDER ALTERING/REMOVING

The program may be used to analyse actual data, by including genotypes at a locus, or to investigate the divergence required for significance in a proposed study. In the latter case assumed base population allele frequencies need to be specified, but no genotypes are required. In the former case base population allele frequencies may either be specified or estimates derived from the data may be used.
