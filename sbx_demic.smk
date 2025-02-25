# -*- mode: Snakemake -*-

import os
import sys
from pathlib import Path


COASSEMBLY_DEMIC_FP = ASSEMBLY_FP / "coassembly_demic"
DEMIC_FP = MAPPING_FP / "demic"

TARGET_DEMIC = DEMIC_FP / "all_PTR.txt"

try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


def get_demic_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_demic":
            return Path(fp)
    raise Error(
        "Filepath for demic not found, are you sure it's installed under extensions/sbx_demic?"
    )


rule all_demic:
    input:
        TARGET_DEMIC,


def zip3l_demic(l1, l2, l3):
    return list(zip(l1, l2, l3))


def coassembly_groups_demic(fp, sample_list):
    if fp == "":
        K = ["all"] * (len(sample_list) * 2)
        V = list(sorted(sample_list)) * 2
        R = [1] * len(sample_list) + [2] * len(sample_list)
        return [K, V, R]
    groups = ruamel.yaml.safe_load(open(str(fp)).read())
    sorted_keys = sorted(groups.keys())
    K = []  # group
    V = []  # sample
    for k in sorted_keys:
        K += [k] * len(groups[k])
        V += groups[k]
    R = [1] * len(V) + [2] * len(V)
    return [K + K, V + V, R]


rule all_coassemble_demic:
    input:
        a=expand(
            COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
            group=list(
                set(
                    coassembly_groups_demic(
                        Cfg["sbx_demic"]["group_file"], Samples.keys()
                    )[0]
                )
            ),
        ),
        b=expand(
            COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_{rp}.fastq",
            zip3l_demic,
            group=coassembly_groups_demic(
                Cfg["sbx_demic"]["group_file"], Samples.keys()
            )[0],
            sample=coassembly_groups_demic(
                Cfg["sbx_demic"]["group_file"], Samples.keys()
            )[1],
            rp=coassembly_groups_demic(Cfg["sbx_demic"]["group_file"], Samples.keys())[
                2
            ],
        ),


rule prep_samples_for_concatenation_paired_demic:
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        r1=temp(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_1.fastq"),
        r2=temp(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_2.fastq"),
    benchmark:
        BENCHMARK_FP / "prep_samples_for_concatenation_paired_demic_{sample}_{group}.tsv"
    log:
        LOG_FP / "prep_samples_for_concatenation_paired_demic_{sample}_{group}.log",
    threads: Cfg["sbx_demic"]["coassembly_threads"]
    conda:
        "envs/sbx_demic_coassembly_env.yml"
    shell:
        """
        pigz -d -p {threads} -c {input.r1} > {output.r1}
        pigz -d -p {threads} -c {input.r2} > {output.r2}
        """


rule combine_groups_paired_demic:
    input:
        rules.all_coassemble_demic.input.b,
    output:
        r1=COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz",
        r2=COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz",
    params:
        w1=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_1.fastq")),
        w2=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_2.fastq")),
    threads: Cfg["sbx_demic"]["coassembly_threads"]
    conda:
        "envs/sbx_demic_coassembly_env.yml"
    resources:
        runtime=120,
    shell:
        """
        cat {params.w1} | pigz -p {threads} > {output.r1}
        cat {params.w2} | pigz -p {threads} > {output.r2}
        """


rule coassemble_paired_demic:
    input:
        r1=COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz",
        r2=COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz",
    output:
        COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
    benchmark:
        BENCHMARK_FP / "coassemble_paired_{group}.tsv"
    log:
        LOG_FP / "coassemble_paired_demic_{group}.log",
    params:
        assembly_dir=str(COASSEMBLY_DEMIC_FP / "{group}"),
    threads: Cfg["sbx_demic"]["coassembly_threads"]
    conda:
        "envs/sbx_demic_coassembly_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -t {threads} -o {params.assembly_dir} 2>&1 | tee {log}
        mv {params.assembly_dir}/final.contigs.fa {output}
        """


rule maxbin:
    input:
        a=expand(
            COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
            group=list(
                set(
                    coassembly_groups_demic(
                        Cfg["sbx_demic"]["group_file"], Samples.keys()
                    )[0]
                )
            ),
        ),
        b=rules.all_coassemble_demic.input.b,
    output:
        directory(COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin"),
    benchmark:
        BENCHMARK_FP / "maxbin.tsv"
    log:
        LOG_FP / "maxbin.log",
    params:
        basename=str(Cfg["all"]["output_fp"]),
        single_genome=str(Cfg["sbx_demic"]["single_genome"]),
        binned_dir=str(COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin"),
        output=str(COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin.001.fasta"),
        maxbin_dir=str(get_demic_path() / "MaxBin_2.2.7_scripts"),
        script=str(get_demic_path() / "MaxBin_2.2.7_scripts" / "run_MaxBin.pl"),
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        find {params.basename}/qc/decontam -iname '*.fastq.gz' > {params.basename}/decontam_list
        mkdir -p {params.binned_dir}

        if [ {params.single_genome} = "True" ] || [ {params.single_genome} = "true" ] || [ {params.single_genome} = "TRUE" ]
        then
            echo "Single genome, skipping MaxBin"
            cp {input.a} {params.output}
            exit 0
        fi

        if command -v MaxBin &> /dev/null
        then
            echo "Using included scripts for MaxBin"
            cd {params.maxbin_dir}
            {params.script} -thread 10 -contig {input.a} \
            -out {params.binned_dir} -reads_list {params.basename}/decontam_list \
            -verbose 2>&1 | tee {log}
        elif command -v run_MaxBin.pl &> /dev/null
        then
            echo "Using built-in scripts for MaxBin"
            run_MaxBin.pl -thread 10 -contig {input.a} \
            -out {params.binned_dir} -reads_list {params.basename}/decontam_list \
            -verbose 2>&1 | tee {log}
        else
            echo "Could not find MaxBin or run_MaxBin.pl in $PATH" > {log}
        fi
        """


rule bowtie2_build:
    input:
        COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin",
    output:
        touch(COASSEMBLY_DEMIC_FP / "max_bin" / ".indexed"),
    params:
        base_dir=str(COASSEMBLY_DEMIC_FP / "max_bin"),
    threads: Cfg["sbx_demic"]["demic_threads"]
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        for f in {params.base_dir}/*.fasta; do
            bowtie2-build --threads {threads} $f $f
        done
        """


rule bowtie2:
    input:
        bin_dir=COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin",
        reads=expand(
            QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
            sample=Samples.keys(),
            rp=Pairs,
        ),
        indexed=COASSEMBLY_DEMIC_FP / "max_bin" / ".indexed",
    output:
        directory(DEMIC_FP / "raw"),
    threads: Cfg["sbx_demic"]["demic_threads"]
    params:
        base_dir=str(COASSEMBLY_DEMIC_FP / "max_bin"),
        reads_dir=str(QC_FP / "decontam"),
    conda:
        "envs/demic_bio_env.yml"
    script:
        "scripts/bowtie2.py"


rule samtools_sort:
    input:
        DEMIC_FP / "raw",
    output:
        directory(DEMIC_FP / "sorted"),
    threads: Cfg["sbx_demic"]["demic_threads"]
    conda:
        "envs/demic_bio_env.yml"
    log:
        str(DEMIC_FP / "logs" / "samtools.error"),
    script:
        "scripts/samtools_sort.py"


rule install_demic:
    output:
        out=DEMIC_FP / ".installed",
    conda:
        "envs/demic_env.yml"
    script:
        "scripts/install_demic.R"


rule run_demic:
    input:
        DEMIC_FP / "sorted",
        COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin",
        DEMIC_FP / ".installed",
    output:
        directory(DEMIC_FP / "DEMIC_OUT"),
    params:
        demic=get_demic_path() / "vendor_demic_v1.0.2" / "DEMIC.pl",
        sam_dir=str(DEMIC_FP / "sorted"),
        fasta_dir=str(COASSEMBLY_DEMIC_FP / "max_bin"),
        keep_all=Cfg["sbx_demic"]["keepall"],
        extras=Cfg["sbx_demic"]["extras"],
    threads: Cfg["sbx_demic"]["demic_threads"]
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "envs/demic_env.yml"
    log:
        LOG_FP / "run_demic.log",
    shell:
        """
        {params.demic} --output_all {params.keep_all} {params.extras} \
        --thread_num {threads} \
        -S {params.sam_dir} -F {params.fasta_dir} \
        -O {output} 2>&1 | tee {log}
        """


rule aggregate_demic:
    input:
        DEMIC_FP / "DEMIC_OUT",
    output:
        DEMIC_FP / "all_PTR.txt",
    shell:
        """
        cat {input}/*.ptr> {output}
        """
