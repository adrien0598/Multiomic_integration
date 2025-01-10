import numpy as np
import pandas as pd
import scanpy as sc

import plotnine as p9

import liana as li
import decoupler as dc
import omnipath as op

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from scipy.io import mmread

#Collectri regulon
net = dc.get_collectri()
PATH_TO_DATA = "/Users/athomas/Documents/PhD_projects/Surgery/Mutliomic_integration/data/"
adata = sc.read_h5ad(PATH_TO_DATA + "GL261_anndata.h5ad")
RNA_counts = mmread(PATH_TO_DATA + 'SCT_counts.mtx').transpose().tocsr().astype('float64')
layers_dict = {'RNA': RNA_counts}
adata.X = layers_dict['RNA']
adata.layers['counts']=adata.X.copy()

sc.pp.filter_cells(adata, min_genes=20)
sc.pp.filter_genes(adata, min_cells=3)

sample_key = 'orig.ident'
groupby = 'cell_types'
condition_key = 'cond'

# PseudoBulk
pdata = dc.get_pseudobulk(
    adata,
    sample_col=sample_key,
    groups_col=groupby,
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000
)

#dc.plot_psbulk_samples(pdata, groupby=[sample_key, groupby], figsize=(11, 4))

## DEA
dea_results = {}
quiet = True

for cell_group in pdata.obs[groupby].unique():
    # Select cell profiles
    ctdata = pdata[pdata.obs[groupby] == cell_group].copy()

    # Obtain genes that pass the edgeR-like thresholds
    # NOTE: QC thresholds might differ between cell types, consider applying them by cell type
    genes = dc.filter_by_expr(ctdata,
                              group=condition_key,
                              min_count=5, # a minimum number of counts in a number of samples
                              min_total_count=10 # a minimum total number of reads across samples
                              )

    # Filter by these genes
    ctdata = ctdata[:, genes].copy()

    # Build DESeq2 object
    # NOTE: this data is actually paired, so one could consider fitting the patient label as a confounder
    dds = DeseqDataSet(
        adata=ctdata,
        design_factors=condition_key,
        ref_level=[condition_key, 'CTL'], # set control as reference
        refit_cooks=True,
        quiet=quiet
    )

    # Compute LFCs
    dds.deseq2()
    # Contrast between stim and ctrl
    stat_res = DeseqStats(dds, contrast=[condition_key, 'RES', 'CTL'], quiet=quiet)
    stat_res.quiet = quiet
    # Compute Wald test
    stat_res.summary()
    # Shrink LFCs
    stat_res.lfc_shrink(coeff='cond_RES_vs_CTL') # {condition_key}_cond_vs_ref

    dea_results[cell_group] = stat_res.results_df

# concat results across cell types
dea_df = pd.concat(dea_results)
dea_df = dea_df.reset_index().rename(columns={'level_0': groupby,'level_1':'index'}).set_index('index')
dea_df.head()

#dea_df[dea_df.isna().any(axis=1)] = float(1)

adata = adata[adata.obs[condition_key]=='RES'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#sc.pl.umap(adata, color=[condition_key, sample_key, groupby], frameon=False, ncols=2)

lr_res = li.mu.df_to_lr(adata,
                           dea_df=dea_df,
                           resource_name='consensus', # NOTE: uses HUMAN gene symbols!
                           expr_prop=float(0.1), # calculated for adata as passed - used to filter interactions
                           groupby="cell_types",
                           stat_keys=['stat', 'pvalue', 'padj'],
                           use_raw=False,
                           complex_col='stat', # NOTE: we use the Wald Stat to deal with complexes
                           verbose=True,
                           return_all_lrs=False,
                           )


lr_res = lr_res.sort_values("interaction_stat", ascending=False, key=abs)
#lr_res.head()
# Let's visualize how this looks like for all interactions  (across all cell types)
lr_res = lr_res.sort_values("interaction_stat", ascending=False)
lr_res['interaction_stat'].hist(bins=50)

df = pd.DataFrame(lr_res)
df.to_csv('lr_res_GL261.csv', index=False)

li.pl.tileplot(liana_res=lr_res,
               fill = 'expr',
               label='padj',
               label_fun = lambda x: '*' if x < 0.05 else np.nan,
               top_n=15,
               orderby = 'interaction_stat',
               orderby_ascending = False,
               orderby_absolute = False,
               source_title='Ligand',
               target_title='Receptor',
               )

plot = li.pl.dotplot(liana_res=lr_res,
                     colour='interaction_stat',
                     size='ligand_pvalue',
                     inverse_size=True,
                     orderby='interaction_stat',
                     orderby_ascending=False,
                     orderby_absolute=True,
                     top_n=10,
                     size_range=(0.5, 4)
                     )

