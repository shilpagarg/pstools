## How to create an initial template for Galaxy tool.

```
$ planemo tool_init --force \
                    --id 'pstools_qv' \
                    --name 'Calculate qv score for prediction (pstools)' \
                    --example_command 'pstools qv seq.fa hic.R1.fastq.gz hic.R2.fastq.gz > qv.txt' \
                    --example_input seq.fa \
                    --example_input hic.R1.fastq.gz \
                    --example_input hic.R2.fastq.gz \
                    --example_output qv.txt \
                    --cite_url 'https://github.com/shilpagarg/pstools' \
                    --help_from_command 'pstools qv --help || :' \
                    --container "ghcr.io/junaruga/pstools:latest"
```
