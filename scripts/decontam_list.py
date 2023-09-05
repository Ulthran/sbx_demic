import ast
import os

if snakemake.params["mapping_fp"]:
    with open(snakemake.output[0], "w") as f_out, open(
        snakemake.params["mapping_fp"]
    ) as f_map, open(snakemake.log[0], "w") as log:
        all_files = [
            fp
            for fp in os.listdir(snakemake.params["decontam_fp"])
            if fp.endswith(".fastq.gz")
        ]

        for line in f_map.readlines():
            if line.split(":")[0] == snakemake.params["group"]:
                sample_names = ast.literal_eval(line.split(":")[1].strip())
                for fp in all_files:
                    if any([n for n in sample_names if n in fp]):
                        log.write(f"Adding {fp}\n")
                        f_out.write(f"{snakemake.params['decontam_fp']}/{fp}\n")
else:
    with open(snakemake.output[0], "w") as f_out, open(snakemake.log[0], "w") as log:
        for fp in [
            fp
            for fp in os.listdir(snakemake.params["decontam_fp"])
            if fp.endswith(".fastq.gz")
        ]:
            log.write(f"Adding {fp}\n")
            f_out.write(f"{snakemake.params['decontam_fp']}/{fp}\n")
