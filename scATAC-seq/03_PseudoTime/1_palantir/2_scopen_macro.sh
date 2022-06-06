#!/usr/bin/env zsh

scopen --input ../save/Macrophages_filtered_peak_bc_matrix  --input-format 10X --output-dir ../save/scOpen --output-prefix Macrophages --max-rho 0.5 --output-format dense --verbose 1 --n-components 30
