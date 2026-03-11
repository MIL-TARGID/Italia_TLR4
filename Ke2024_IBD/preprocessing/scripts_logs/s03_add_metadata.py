import scanpy as sc
import pandas as pd
import os
import argparse
import datetime

# env is norm

### logging
def log(msg, log_file):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    full_msg = f"{timestamp} {msg}"
    print(full_msg)
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(full_msg + '\n')

### add metadata
def add_metadata(input_directory_path, input_metadata_path):

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(
        os.getcwd(),  # logging where called
        f"add_metadata_log_{timestamp}.txt"
    )
    log("🚀 Starting metadata update!\n", log_file)

    # load metadata csv
    metadata_df = pd.read_csv(input_metadata_path)
    log(f"📄 Loaded metadata file: {input_metadata_path}\n", log_file)
    log("🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍🤍\n", log_file)

    # loop over h5ad files
    for file in os.listdir(input_directory_path):
        if file.endswith(".h5ad"):
            sample_path = os.path.join(input_directory_path, file)
            sample_name = os.path.splitext(file)[0]  # filename without .h5ad
            log(f"🤖 Processing file: {file}", log_file)

            adata = sc.read_h5ad(sample_path)

            # add 'sample' column to obs
            adata.obs['sample'] = sample_name
            log(f"✅ Added new column 'sample' with '{sample_name}'", log_file)

            # prepend obs_names with sample_name
            adata.obs_names = [f"{sample_name}_{cell}" for cell in adata.obs_names]
            log(f"✅ Modified obs_names barcodes with sample name prefix", log_file)

            # match metadata rows
            matched = metadata_df[metadata_df['Sample'] == sample_name]

            matched_row = matched.iloc[0]

            # automatic metadata column detection
            metadata_fields = [col for col in metadata_df.columns if col != "Sample"]

            for col in metadata_fields:
                adata.obs[col] = matched_row[col]
                log(f"✅ Added metadata column: '{col}'", log_file)

            adata.write(sample_path)
            log(f"💾 Saved updated file: {file}", log_file)
            log(f"📊 AnnData summary:\n{adata}", log_file)
            log(f"🔍 Sample obs preview:\n{adata.obs.head(2)}", log_file)
            log(f"🔍 Sample var preview:\n{adata.var.head(2)}\n", log_file)
            log("❤️🧡💛💚💙💜❤️🧡💛💚💙💜❤️\n", log_file)

    log("🎉 All metadata updated! \n", log_file)
    log(f"💾 Log saved to: {log_file}\n", log_file)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add metadata to anndata files using .csv file")
    parser.add_argument("input_directory", type=str, help="Path to directory with .h5ad files")
    parser.add_argument("input_metadata", type=str, help="Path to metadata .csv file")

    args = parser.parse_args()
    add_metadata(args.input_directory, args.input_metadata)