# module load $scvi
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad

# Options
anndataf = "/lustre/scratch126/humgen/teams_v2/anderson/cc59/data/scRNA/adata_celltypist_annotated_filtered_3.h5ad"
method = "sum"
cell_annot_col = "Category"
ind_col = "Batch123456_Genotyping_ID"
pseudobulk_outf = f"/lustre/scratch126/humgen/teams_v2/anderson/cc59/data/scRNA/pseudobulk_adata_{cell_annot_col}-by-{ind_col}-d{method}.h5ad"

# Function
def pseudobulk_by_label(adata, groupby='pool_participant', label='Celltypist:IBDverse_eqtl:predicted_labels',
                        min_cells=5, min_samples=30, layer='log1p_cp10k', method='mean',
                        min_samples_per_gene=0.0, nhvgs=None):
    """Pseudobulk expression by label (cell type) and sample.
    Returns a dict keyed by label value, values are DataFrames (samples x genes).
    Only includes labels with at least min_samples samples.
    """
    print(f"..Fetching layer '{layer}'")
    if layer in adata.layers:
        expr = adata.layers[layer]
    else:
        raise ValueError(f"Layer '{layer}' not found in adata.layers")
    print("..Copying obs")
    obs = adata.obs.copy()
    obs['_cell_idx'] = np.arange(len(obs))
    labels = obs[label].dropna().unique()
    print(f"..Found {len(labels)} labels in '{label}'")
    pseudobulk_dict = {}
    for lab in labels:
        print(f"..Processing label: {lab}")
        mask = (obs[label] == lab).values
        obs_sub = obs[mask].copy()
        cell_counts = obs_sub.groupby(groupby).size()
        valid_samples = cell_counts[cell_counts >= min_cells].index
        if len(valid_samples) < min_samples:
            print(f"....Skipping '{lab}': only {len(valid_samples)} samples pass "
                  f"min_cells={min_cells} (need {min_samples})")
            continue
        obs_sub = obs_sub[obs_sub[groupby].isin(valid_samples)]
        # Only densify the cells needed for this label, not the whole matrix
        print(f"....Densifying {len(obs_sub)} cells x {adata.shape[1]} genes")
        expr_sub = expr[obs_sub['_cell_idx'].values, :]
        if hasattr(expr_sub, 'toarray'):
            expr_sub = expr_sub.toarray()
        else:
            expr_sub = np.asarray(expr_sub)
        sample_ids = obs_sub[groupby].values
        unique_samples = np.unique(sample_ids)
        pseudobulk_expr = np.zeros((len(unique_samples), expr_sub.shape[1]))
        for i, samp in enumerate(unique_samples):
            samp_mask = sample_ids == samp
            if method == 'mean':
                pseudobulk_expr[i, :] = expr_sub[samp_mask, :].mean(axis=0)
            elif method == 'sum':
                pseudobulk_expr[i, :] = expr_sub[samp_mask, :].sum(axis=0)
            else:
                raise ValueError(f"Unknown method: {method}")
        pseudobulk_df = pd.DataFrame(pseudobulk_expr, index=unique_samples, columns=adata.var_names)
        if min_samples_per_gene > 0:
            n_samples = len(unique_samples)
            min_sample_count = int(np.ceil(min_samples_per_gene * n_samples))
            genes_expressed = (pseudobulk_df > 0).sum(axis=0) >= min_sample_count
            pseudobulk_df = pseudobulk_df.loc[:, genes_expressed]
            print(f"......Kept {genes_expressed.sum()}/{len(genes_expressed)} genes expressed in "
                  f">={min_samples_per_gene:.1%} of samples")
        if nhvgs is not None:
            gene_vars = pseudobulk_df.var(axis=0)
            if isinstance(nhvgs, float) and 0 < nhvgs < 1:
                n_to_select = int(np.ceil(nhvgs * len(gene_vars)))
            elif isinstance(nhvgs, int) and nhvgs > 0:
                n_to_select = min(nhvgs, len(gene_vars))
            else:
                raise ValueError(f"nhvgs must be a positive int or float between 0 and 1, got {nhvgs}")
            top_var_genes = gene_vars.nlargest(n_to_select).index
            pseudobulk_df = pseudobulk_df[top_var_genes]
            print(f"......Selected {n_to_select}/{len(gene_vars)} highly variable genes by variance")
        pseudobulk_dict[lab] = pseudobulk_df
        print(f"....Done '{lab}': {pseudobulk_df.shape[0]} samples x {pseudobulk_df.shape[1]} genes")
    return pseudobulk_dict

def combine_pseudobulk(pseudobulk_dict, ind_col, cell_annot_col):
    """Vertically combine a {label: DataFrame(samples x genes)} dict into one DataFrame.
    Rownames become "<sample>-<label>". Genes that aren't shared across every
    label are kept (missing entries filled with NaN via outer join on columns).
    Also returns a rowname -> (sample, label) map as a DataFrame (columns
    `ind_col`, `cell_annot_col`) built from the dict key/index directly, so
    hyphens inside sample IDs or labels can't corrupt a later re-parse of the
    combined rowname.
    """
    print(f"..Combining {len(pseudobulk_dict)} labels")
    combined_parts = []
    meta_parts = []
    for lab, df in pseudobulk_dict.items():
        samples = df.index
        new_index = [f"{sample}-{lab}" for sample in samples]
        df_relabelled = df.copy()
        df_relabelled.index = new_index
        combined_parts.append(df_relabelled)
        meta_parts.append(pd.DataFrame({ind_col: samples, cell_annot_col: lab}, index=new_index))
        print(f"....Added '{lab}': {df_relabelled.shape[0]} rows")
    combined = pd.concat(combined_parts, axis=0, join='outer')
    meta = pd.concat(meta_parts, axis=0)
    print(f"..Combined: {combined.shape[0]} rows x {combined.shape[1]} genes")
    return combined, meta

########
# Load anndata and pseudobulk
########
print(f"Loading h5ad: {anndataf}")
adata = sc.read_h5ad(anndataf)
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

pb = pseudobulk_by_label(
            adata,
            groupby=ind_col,
            label=cell_annot_col,
            min_cells=5,
            min_samples=30,
            layer="counts",
            method=method,
            min_samples_per_gene=0,
            nhvgs=None, # Taking all genes for now, let the pipeline handle filtering
        )

print(f"Finished pseudobulking: {len(pb)}/{adata.obs[cell_annot_col].nunique()} labels retained")

pb_combined, pb_meta = combine_pseudobulk(pb, ind_col, cell_annot_col)


########
# Manually create some obs
########

# Demographic covs
cols_to_keep = ["sex", "age"]
demo_covs = adata.obs[[ind_col, cell_annot_col] + cols_to_keep].copy().reset_index(drop=True).drop_duplicates(subset=[ind_col, cell_annot_col])
demo_covs = demo_covs.set_index(demo_covs[ind_col].astype(str) + "-" + demo_covs[cell_annot_col].astype(str))
pb_meta = pb_meta.join(demo_covs[cols_to_keep], how="left")


########
# Manually create the anndata
########
adata_pb = ad.AnnData(
    X=pb_combined.values,
    obs=pb_meta,
    var=adata.var[['gene_symbols', 'feature_types', 'gene_group__mito_transcript',
       'gene_group__mito_protein', 'gene_group__ribo_protein',
       'gene_group__ribo_rna', 'mt']] # Keep these from the original var
)


#########
# Re-derive some technical covariates at the pseubulk level
#########
print("..Computing pseudobulk QC metrics (mito %, genes expressed)")
sc.pp.calculate_qc_metrics(adata_pb, qc_vars=['gene_group__mito_transcript'], percent_top=None, log1p=False, inplace=True)
# adds to adata_pb.obs: n_genes_by_counts, total_counts, total_counts_mt, pct_counts_mt

# Save
print(f"Saving pseudobulk anndata to: {pseudobulk_outf}")
adata_pb.write_h5ad(pseudobulk_outf)