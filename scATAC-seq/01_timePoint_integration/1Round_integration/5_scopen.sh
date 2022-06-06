#!/bin/bash
sample=$1
scopen --input ../../1Round_Peakcalling/data/${sample}/filtered_peak_bc_matrix --input-format 10X --output-dir ./ --output-prefix Heart --max-rho 0.5 --output-format dense --verbose 1 --n-components 30

