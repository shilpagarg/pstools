# pstools: Toolkit for full phased sequences of diploid genomes 

Routine production of phased sequences (haplotypes) of genomes are important to study variation occurring in complex cancer and healthy samples. A combination of long accurate reads (HiFi) and chromosome-scale sequencing (Hi-C) technologies are enabling the routine production of fully phased sequences. However, standard computational approaches have limited capability to produce phased sequences that are end-to-end chromosomes in less than a day without any collapse in haplotypes. Here, we propose a fast and accurate graph-based method (pstools) that jointly leverage phasing and sequencing information from HiFi and Hi-C to produce high-quality phased sequences. We benchmarked this method on the normal sample (HG002) to produce sequence continuity N50 >130 Mb, base quality >Q50, and completeness of 3Gb of each haplotype, by comparing against GIAB ground truth. 


