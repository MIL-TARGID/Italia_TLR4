import scanpy as sc
import pandas as pd
import os
import argparse
import datetime
import logging
import numpy as np

import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

from rpy2.robjects import globalenv
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from anndata2ri import py2rpy
from anndata2ri import converter as ann_converter

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

import rpy2.rinterface_lib.callbacks
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.INFO)
rpy2.rinterface_lib.callbacks.consolewrite_print = lambda x: print(x, end="")
rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda x: print(x, end="")


# R code using rpy2
r = ro.r

# syntax for R magic
r('''
library(SoupX)
library(SingleCellExperiment)
''')

### logging helper
def log(msg, log_file):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    full_msg = f"{timestamp} {msg}"
    print(full_msg)
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(full_msg + '\n')


def ambient(input_path):

    # output path
    parent_path = os.path.dirname(input_path)
    output_path = os.path.join(parent_path, "out_ambientRNA") # moves path one up from input_path
    os.makedirs(output_path, exist_ok=True) # create if not existing

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(
        os.getcwd(),  # logging where called
        f"ambient_removal_log_{timestamp}.txt"
    )
    log("🚀 Starting ambient RNA removal!\n", log_file)

    log(f"📂 Output path: {output_path}", log_file)

    # check/create output subdirs
    output_metrics_path = os.path.join(output_path, "metrics")
    os.makedirs(output_metrics_path, exist_ok=True)
    log(f"💾 SoupX metrics and plots will be saved here: {output_metrics_path}", log_file)

    output_anndata_path = os.path.join(output_path, "anndata")
    os.makedirs(output_anndata_path, exist_ok=True)
    log(f"💾 SoupX-corrected AnnData files will be saved here: {output_anndata_path}", log_file)

    # list all subdirs (samples)
    directories = [f for f in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, f))]
    log(f"📂 Found subdirectories (samples): {directories}\n", log_file)
    log("🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍\n", log_file)

    # iterate over subdirs
    for direc in directories:
        direc_path = os.path.join(input_path, direc)
        sample = direc
        df_path = os.path.join(direc_path, "outs", "analysis", "clustering", "gene_expression_graphclust", "clusters.csv")
        adata_path = os.path.join(direc_path, "outs", "filtered_feature_bc_matrix.h5")
        adata_raw_path = os.path.join(direc_path, "outs", "raw_feature_bc_matrix.h5")

        # check files exist
        if not all(os.path.exists(p) for p in [df_path, adata_path, adata_raw_path]):
            log(f"⚠️ Missing one or more files for sample {direc}", log_file)
            continue
        
        log(f"🤖 Processing sample: {sample}\n", log_file)

        # load cellranger clustering csv
        df = pd.read_csv(df_path)
        # series: index = Barcode, values = Cluster, type = category
        soupx_groups = pd.Series(df["Cluster"].astype(str).values, index=df["Barcode"]).astype("category")
        del(df)

        # callranger filtered h5
        adata = sc.read_10x_h5(filename=adata_path)

        # duplicates in var names
        old_names = adata.var_names.copy()
        duplicated_vars = old_names[old_names.duplicated()]

        if len(duplicated_vars) > 0:
            log(f"⚠️ Found {len(duplicated_vars)} duplicated var names:", log_file)
            log(f"😲 Duplicates: {duplicated_vars.tolist()}", log_file)
        else:
            log(f"✅ No duplicated var names found", log_file)

        # make var unique again
        adata.var_names_make_unique()

        # which var changed
        new_names = adata.var_names
        changed = [(old, new) for old, new in zip(old_names, new_names) if old != new]

        if changed:
            log(f"🧩 Renamed duplicated vars:", log_file)
            for old, new in changed:
                log(f"   {old} → {new}", log_file)
            log(f"✅ Total renamed vars: {len(changed)}\n", log_file)
        else:
            log(f"✅ No var names were changed", log_file)

        # extract cells(barcodes), genes and X transposed
        cells = adata.obs_names
        genes = adata.var_names

        # load cellranger raw h5
        adata_raw = sc.read_10x_h5(filename=adata_raw_path)
        adata_raw.var_names_make_unique()

        with localconverter(ro.default_converter + pandas2ri.converter):
            globalenv["sample"] = sample
            globalenv["output_metrics_path"] = output_metrics_path
            globalenv["soupx_groups"] = ro.conversion.py2rpy(soupx_groups)
            globalenv["genes"] = ro.conversion.py2rpy(genes)
            globalenv["cells"] = ro.conversion.py2rpy(cells)
 
        # pass full anndata to R
        with localconverter(ro.default_converter + ann_converter):
            globalenv["adata"] = adata
            globalenv["adata_raw"] = adata_raw


        r('''

          # row and column names of data
          rownames(adata) = genes
          colnames(adata) = cells

          # table of counts and table of droplets
          data <- assay(adata, "X")
          data_raw <- assay(adata_raw, "X")

          # generate soupchannel object for soupx
          sc = SoupChannel(data_raw, data)

          # add metadata to the object
          soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
          sc = setSoupProfile(sc, soupProf)

          # set cluster information in obj
          sc = setClusters(sc, soupx_groups)

          # estimate contamination fraction + plot Rho auto estimates
          png(paste0(output_metrics_path, "/", sample, "_soupx_autoEstCont_plot.png"), width=1000, height=800, res=200)
          sc <- autoEstCont(sc, doPlot=TRUE)
          dev.off()

          # plot marker distribution
          png(paste0(output_metrics_path, "/", sample, "_soupx_plotMarkerDistribution.png"), width=1600, height=1200, res=200)
          markdis <- plotMarkerDistribution(sc)
          print(markdis)
          dev.off()

          # save metadata and contamination fraction
          metadata <- sc$metaData
          estrho <- sc$fit$rhoEst
          fit_dd <- sc$fit$dd 
          markers_used <- sc$fit$markersUsed
          soupProfile <- sc$soupProfile

          # save this value as plain text - no header, no row names, no quotes
          write.table(
            estrho,
            file = paste0(output_metrics_path, "/", sample, "_soupx_rho.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
          )

          write.csv(
            metadata,
            file = paste0(output_metrics_path, "/", sample, "_soupx_metadata.csv"),
            row.names = TRUE,
            col.names = TRUE
          )

          write.csv(
            fit_dd,
            file = paste0(output_metrics_path, "/", sample, "_soupx_fit_dd.csv"),
            row.names = TRUE,
            col.names = TRUE
          )

          write.csv(
            markers_used,
            file = paste0(output_metrics_path, "/", sample, "_soupx_markers_used.csv"),
            row.names = TRUE,
            col.names = TRUE
          )

          write.csv(
            soupProfile,
            file = paste0(output_metrics_path, "/", sample, "_soupx_soupProfile.csv"),
            row.names = TRUE,
            col.names = TRUE
          )

          # corrected table of counts - round to integer so downstream tools dont cry + less memory
          out = adjustCounts(sc, roundToInt = TRUE)
          ''')

        with localconverter(ro.default_converter + ann_converter):
            out_corrected = ro.conversion.rpy2py(ro.globalenv["out"])

        adata.layers["counts"] = adata.X
        adata.layers["soupX_counts"] = out_corrected.T
        adata.X = adata.layers["soupX_counts"]
        adata.write(f"{output_anndata_path}/{sample}.h5ad")
        log(f"💾 {sample} AnnData with 'counts' and 'soupX_counts' layers saved: {output_anndata_path}/{sample}.h5ad", log_file)
        log(f"📊 {sample} Anndata summary:\n{adata}", log_file)

        # sum across cells axis=0 - total gene expression
        original = np.asarray(adata.layers["counts"].sum(axis=0)).ravel()
        corrected = np.asarray(adata.layers["soupX_counts"].sum(axis=0)).ravel()
        diff = original - corrected
        gene_diff = pd.DataFrame({"gene": adata.var_names, "removed_counts": diff})
        gene_diff = gene_diff.sort_values(by="removed_counts", ascending=False)
        gene_diff.to_csv(f'{output_metrics_path}/{sample}_genes_removed_counts.csv', index=False)
        log(f"📊 {sample} table of removed gene counts saved: {output_metrics_path}/{sample}_genes_removed_counts.csv", log_file)

        # sum across counts per cell axis=1
        original = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
        corrected = np.asarray(adata.layers["soupX_counts"].sum(axis=1)).ravel()
        diff = original - corrected
        gene_diff = pd.DataFrame({"gene": adata.obs_names, "removed_counts": diff})
        gene_diff = gene_diff.sort_values(by="removed_counts", ascending=False)
        gene_diff.to_csv(f'{output_metrics_path}/{sample}_counts_per_cell_removed.csv', index=False)
        log(f"📊 {sample} table of removed counts per cell saved: {output_metrics_path}/{sample}_counts_per_cell_removed.csv", log_file)
        log(f"✅ Finished processing sample: {sample}\n", log_file)
        log("❤️🧡💛💚💙💜❤️🧡💛💚💙💜❤️\n", log_file)

    log("🎉 All ambientRNA removal completed!", log_file)
    log(f"💾 Log saved to: {log_file}\n", log_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Remove ambient RNA with SoupX from CellRanger out directory")
    parser.add_argument("input", type=str, help="Path to master CellRanger output directory with sample subdirectories")

    args = parser.parse_args()
    ambient(args.input)