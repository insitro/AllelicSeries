## Version 1.0.0
* Added ability to run COAST starting from summary statistics:
	- `ASBTSS` runs the allelic series burden test from summary statistics.
	- `ASKATSS` runs the allelic series SKAT test from summary statistics.
	` `COASTSS` runs the omnibus coding-variant allelic series test from summary statistics. 
* Added a `CalcSumstats` function to generate summary statistics from individual-level data. 

## Version 0.6.0
* Replaced the `CountAlleles` function with a similar (but faster) `Counts` function that returns counts the total number of alleles, variants, and carriers by variant class.
* Updated the formatting of the results output by the main `COAST` function.

## Version 0.5.0

* Added option (`min_mac`) to filter the variant set to only include those variants having at least a minimum minor allele count (10 is recommended).
* Added a function (`CountAlleles`) to count the number of alleles of each variant category present in the genotype matrix. Also allows for counting the number of carriers of each type of allele.
* By default, `COAST` now reports the number of alleles of each variant category that contributed to the test. 

