# pstools: Toolkit for fully phased sequences

Routine production of phased sequences (haplotypes) of genomes are important to study variation occurring in complex cancer and healthy samples. We developed a novel graph-based algorithm to integrate HiFi and Hi-C data types to produce haplotypes at the base-level resolution for routine clinical research. When benchmarking our method on healthy human genomes, we produced significantly high-quality genomes with a sequence continuity NG50 >130 Mb, switch/hamming error rates <1.5% and a completeness of >6.0 Gb, and an order of magnitude faster process (only <12 hours). This simple approach will facilitate improvement in understanding of the mechanisms of clinical diseases.

## Installation
```sh
(requiring g++ and zlib)
git clone https://github.com/shilpagarg/pstools.git
cd pstools && make
```

## Execution
Once you compiled the repo, you will see binary file `pstools` in the same directory.
Now you are ready to go to produce fully phased sequences.

```sh
# Use Hi-C data and node sequences of hifiasm graph (awk '/^S/{print ">"$2;print $3}' hifiasm_r_utg.gfa > hifiasm_r_utg.fa)
pstools hic_mapping -t32 -o <map.out> <hifiasm_r_utg.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>

# Use Hi-C mapped reads through hifiasm graph to produce fully phased sequences
pstools resolve_haplotypes -t32 -i true <map.out> <hifiasm_r_utg.gfa> <out>
# where `map.out` is the file from above process, `hifiasm_r_utg.gfa` is the hifiasm r_utg graph and `out` is the output directory name.
This command produces fully phased sequences in `pred_hap1.fa` and `pred_hap2.fa`.  
```

### Results
Table shows the benchmarking results of HG002 using OmniC or Arima genomics data: 

|<sub>Hi-C dataset<sub>|<sub>Size<sub>|<sub>CPU time<sub>|<sub> NG50<sub>|<sub> QV<sub>|
|:---------------|-----:|--------------------:|:----------|-------:|
|<sub>[HG002]</sub>|<sub>~6.1Gb</sub> |<sub>~5h</sub> |<sub>~132 Mb</sub>|<sub>~Q50</sub>|


[HG0002-data]: http://dovetail-omnic.s3-website-us-west-2.amazonaws.com/ or https://www.biorxiv.org/content/10.1101/810341v1, and the resultant phased sequences are available at: s3://pstools/
 
 ```sh
 For further experiments, please run `experiments.sh`.
```
 
### Limitations
1. At the current stage, the phased sequences don't contain centromeric regions.
2. The UL nanopore data is not used.
3. Need to be tested on trio-hifiasm graph
 
Note that pstools is under active development. For any issues or suggestions, please create an issue or a pull request. 

