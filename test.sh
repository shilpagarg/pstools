#!/bin/bash

wget https://test1-sg.s3.us-east-2.amazonaws.com/assembly1.gfa
wget https://test1-sg.s3.us-east-2.amazonaws.com/hic_name_connection1.output

wget https://github.com/shilpagarg/pstools/releases/download/v0.1/pstools
chmod 777 pstools

./pstools resolve_haplotypes -t32 hic_name_connection1.output assembly1.gfa ./
