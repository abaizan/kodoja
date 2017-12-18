import os

old_path = os.path.join(os.environ["HOME"], "Scratch/kraken_taxonomy/")
wanted = set(["137758", "946046", "12227"])  # taxids for our three viruses

for filename in ["nucl_est.accession2taxid",
                 "nucl_gb.accession2taxid",
                 "nucl_gss.accession2taxid",
                 "nucl_wgs.accession2taxid"]:
    with open(os.path.join(old_path, filename)) as input_handle:
        with open(filename, "w") as output_handle:
            for line in input_handle:
                taxid = line.split("\t")[2]
                if taxid in wanted:
                    output_handle.write(line)
    print(filename + " done")
