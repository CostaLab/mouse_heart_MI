#!/usr/bin/env python
# coding: utf-8
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import pickle
import os
import matplotlib
import matplotlib.pyplot as plt
import palantir
import harmony


topcelltype="Fibroblasts"


## not use will lead 0 dimension for ms_data
dr = pd.read_csv(f"../save/scOpen/{topcelltype}_barcodes.txt", sep="\t", index_col=0).T
data = sc.read_csv(f"../save/scOpen/{topcelltype}.txt", delimiter="\t").transpose()
refadata = sc.read_h5ad("../save/h5ad")
data.obs['celltype'] = refadata.obs['celltype'][data.obs_names]
sc.pp.normalize_per_cell(data)
sc.pp.highly_variable_genes(data, n_top_genes=2000, flavor='cell_ranger')
pca_projections, _ = palantir.utils.run_pca(data, use_hvg=False)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
fdl = harmony.plot.force_directed_layout(dm_res['kernel'], data.obs_names)

nms_fib = data.obs['celltype'][data.obs['celltype'] == "Fibroblasts"].index
x = fdl.loc[nms_fib][['y']].idxmin()
start_cell = x['y']

pdsc = pd.Series(["othercells"]*len(data.obs['time']), index =data.obs['time'].index)
data.obs['root_cell'] = pdsc
data.obs['root_cell'][start_cell] = "rootcell"
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500)
data.obs['time'] = refadata.obs['time.ident'][data.obs_names]
data.obs['gender'] = refadata.obs['gender.ident'][data.obs_names]


obj_dict = {
    "name": "Fibroblasts",
    "pca_projections": pca_projections,
    "data": data,
    "dm_res": dm_res,
    "ms_data": ms_data,
    "fdl": fdl,
    "pr_res": pr_res,
    "start_cell": start_cell
}
pickle.dump(obj_dict, file=open("../save/Palantir_fib.pickle", 'wb'))

