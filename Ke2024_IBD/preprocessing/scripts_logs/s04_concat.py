import anndata as ad
import scanpy as sc
import pandas as pd
import os
from datetime import datetime

def log(log_file, message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    full_message = f"{timestamp} {message}\n"
    print(full_message.strip())
    log_file.write(full_message)
    log_file.flush()

def concatenate_anndata(input_path):

    dataset_dir = input_path
    sample_names = [f for f in os.listdir(dataset_dir) if f.endswith('.h5ad')]
    adata_list = []

    # create/open log file, append

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(
        os.getcwd(),  # logging where called
        f"concatenate_log_{timestamp}.txt"
    )

    with open(log_path, "a") as log_file:
        log(log_file, f"🚀 Concatenation in directory: {dataset_dir}\n")

        for sample_name in sample_names:
            sample = os.path.join(dataset_dir, sample_name)
            single_adata = sc.read_h5ad(sample)

            # check for duplicates in var names
            old_names = single_adata.var_names.copy()
            duplicated_vars = old_names[old_names.duplicated()]

            if len(duplicated_vars) > 0:
                log(log_file, f"⚠️ Found {len(duplicated_vars)} duplicated var names in {sample_name}:")
                log(log_file, f"😲 Duplicates: {duplicated_vars.tolist()}")
            else:
                log(log_file, f"✅ No duplicated var names found in {sample_name}")

            # make var unique again
            single_adata.var_names_make_unique()

            # which var changed
            new_names = single_adata.var_names
            changed = [(old, new) for old, new in zip(old_names, new_names) if old != new]

            if changed:
                log(log_file, f"🧩 Renamed duplicated vars in {sample_name}:")
                for old, new in changed:
                    log(log_file, f"   {old} → {new}")
                log(log_file, f"✅ Total renamed vars: {len(changed)}\n")
            else:
                log(log_file, f"✅ No var names were changed in {sample_name}")


            adata_list.append(single_adata)
            log(log_file, f"🧩 Added {sample_name} to the AnnData list!")

        if not adata_list:
            log(log_file, "❌ No .h5ad files found!")
            return

        log(log_file, f"⚙️ Concatenating {len(adata_list)}\n")

        #concatenate the list of .h5ad files
        adata = ad.concat(adata_list, join="outer")

        ### because we are loosing vars add them separately

        # grab all var dfs from our dictionary
        all_var = [x.var for x in adata_list]
        # concatenate them, keep all
        all_var = pd.concat(all_var, join="outer")
        # remove duplicates by index (keep first occurrence)
        all_var = all_var[~all_var.index.duplicated(keep="first")]
        # add var to adata make sure the order is in accordance
        adata.var = all_var.loc[adata.var_names]

        # concatenated anndata will be named as the master directory, here anndata.h5ad
        dataset_dir = input_path.rstrip("/")  # remove trailing slash if present
        out_file = os.path.join(
            dataset_dir,
            f"{os.path.basename(dataset_dir)}.h5ad"
        )

        adata.write(out_file)

        log(log_file, f"💾 Saved concatenated AnnData to {out_file}")
        log(log_file, f"📊 Final AnnData summary:\n{adata}")

        log(log_file, "❤️🧡💛💚💙💜❤️🧡💛💚💙💜❤️\n")

        log(log_file, f"🔍 Sample obs preview:\n{adata.obs.head(2)}\n")
        log(log_file, f"🔍 Sample var preview:\n{adata.var.head(2)}\n")

        log(log_file, "❤️🧡💛💚💙💜❤️🧡💛💚💙💜❤️\n")

        log(log_file, "🎉 Concatenation completed!")
        log(log_file, f"💾 Log saved to: {log_path}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Concatenating individual samples h5ad files")
    parser.add_argument("input", type=str, help="Path to the directory containing h5ad files")

    args = parser.parse_args()
    concatenate_anndata(args.input)