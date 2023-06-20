# -*- mode: Snakemake -*-

import os
import ruamel.yaml
import sys
from pathlib import Path

COASSEMBLY_DEMIC_FP = ASSEMBLY_FP / "coassembly_demic"
DEMIC_FP = Path(Cfg["all"]["output_fp"]) / "demic"
TARGET_DEMIC = [str(MAPPING_FP / "demic" / "DEMIC_OUT" / "all_PTR.txt")]


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


def get_demic_path() -> Path:
    demic_path = Path(sunbeam_dir) / "extensions" / "sbx_demic"
    if demic_path.exists():
        return demic_path
    raise Error(
        "Filepath for demic not found, are you sure it's installed under extensions/sbx_demic?"
    )


rule all_demic:
    input:
        TARGET_DEMIC,


def zip3l(l1, l2, l3):
    return list(zip(l1, l2, l3))


def coassembly_groups(fp, sample_list):
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
            str(COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa"),
            group=list(
                set(
                    coassembly_groups(
                        Cfg["coassembly_demic"]["group_file"], Samples.keys()
                    )[0]
                )
            ),
        ),
        b=expand(
            str(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_{rp}.fastq"),
            zip3l,
            group=coassembly_groups(
                Cfg["coassembly_demic"]["group_file"], Samples.keys()
            )[0],
            sample=coassembly_groups(
                Cfg["coassembly_demic"]["group_file"], Samples.keys()
            )[1],
            rp=coassembly_groups(
                Cfg["coassembly_demic"]["group_file"], Samples.keys()
            )[2],
        ),


rule prep_samples_for_concatenation_paired_demic:
    input:
        r1=str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
        r2=str(QC_FP / "decontam" / "{sample}_2.fastq.gz"),
    output:
        r1=temp(str(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_1.fastq")),
        r2=temp(str(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_2.fastq")),
    benchmark:
        BENCHMARK_FP / "prep_samples_for_concatenation_paired_{sample}_{group}.tsv"
    log:
        LOG_FP / "prep_samples_for_concatenation_paired_demic_{sample}_{group}.log",
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "sbx_coassembly_env.yml"
    shell:
        """
        pigz -d -p {threads} -c {input.r1} > {output.r1}
        pigz -d -p {threads} -c {input.r2} > {output.r2}
        """


rule combine_groups_paired_demic:
    input:
        rules.all_coassemble_demic.input.b,
    output:
        r1=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz"),
        r2=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz"),
    params:
        w1=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_1.fastq")),
        w2=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_2.fastq")),
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "sbx_coassembly_env.yml"
    resources:
        runtime=120,
    shell:
        """
        cat {params.w1} | pigz -p {threads} > {output.r1}
        cat {params.w2} | pigz -p {threads} > {output.r2}
        """


rule coassemble_paired_demic:
    input:
        r1=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz"),
        r2=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz"),
    output:
        str(COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa"),
    benchmark:
        BENCHMARK_FP / "coassemble_paired_{group}.tsv"
    log:
        LOG_FP / "coassemble_paired_demic_{group}.log",
    params:
        assembly_dir=str(COASSEMBLY_DEMIC_FP / "{group}"),
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "sbx_coassembly_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -t {threads} -o {params.assembly_dir} 2>&1 | tee {log}
        mv {params.assembly_dir}/final.contigs.fa {output}
        """


rule decontam_list:
    input:
        expand(
            QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
            sample=Samples.keys(),
            rp=Pairs,
        ),
    output:
        DEMIC_FP / "decontam_list_{group}.txt",
    benchmark:
        BENCHMARK_FP / "decontam_list_{group}.tsv"
    log:
        LOG_FP / "decontam_list_{group}.log",
    params:
        decontam_fp=QC_FP / "decontam",
        mapping_fp=Cfg["coassembly_demic"]["group_file"],
        group="{group}",
    script:
        "scripts/decontam_list.py"


rule maxbin:
    input:
        a=COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
        b=rules.all_coassemble_demic.input.b,
        decontam_list=DEMIC_FP / "decontam_list_{group}.txt",
    output:
        directory(DEMIC_FP / "maxbin" / "{group}/"),
    benchmark:
        BENCHMARK_FP / "maxbin_{group}.tsv"
    log:
        LOG_FP / "maxbin_{group}.log",
    params:
        maxbin_dir=str(get_demic_path() / "MaxBin_2.2.7_scripts"),
        script=str(get_demic_path() / "MaxBin_2.2.7_scripts" / "run_MaxBin.pl"),
        out_dir=str(DEMIC_FP / "maxbin" / "{group}"),
    conda:
        "sbx_demic_bio_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        mkdir {output}
        cd {params.maxbin_dir}
        
        if command -v MaxBin &> /dev/null
        then
            {params.script} -thread 10 -contig {input.a} \
            -out {params.out_dir} -reads_list {input.decontam_list} \
            -verbose 2>&1 | tee {log}
        elif command -v run_MaxBin.pl &> /dev/null
        then
            run_MaxBin.pl -thread 10 -contig {input.a} \
            -out {params.out_dir} -reads_list {input.decontam_list} \
            -verbose 2>&1 | tee {log}
        else
            echo "Could not find MaxBin or run_MaxBin.pl in $PATH"
        fi
        """


rule bowtie2_build:
    input:
        COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
    output:
        [
            DEMIC_FP / "bowtie2" / "{group}" / ("contigs" + ext)
            for ext in [
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ]
        ],
    benchmark:
        BENCHMARK_FP / "bowtie2_build_{group}.tsv"
    log:
        LOG_FP / "bowtie2-build_{group}.log",
    params:
        basename=str(DEMIC_FP / "bowtie2" / "{group}" / "contigs"),
    threads: 4
    conda:
        "sbx_demic_bio_env.yml"
    resources:
        mem_mb=8000,
        runtime=240,
    shell:
        """
        mkdir -p {params.basename}
        bowtie2-build --threads {threads} {input} {params.basename} 2>&1 | tee {log}
        """


# Run bowtie2 with index
rule bowtie2:
    input:
        DEMIC_FP / "bowtie2" / "{group}" / "contigs.1.bt2",
        rp1=expand(
            str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
            sample=Samples.keys(),
        ),
        rp2=expand(
            str(QC_FP / "decontam" / "{sample}_2.fastq.gz"),
            sample=Samples.keys(),
        ),
    output:
        temp(MAPPING_FP / "demic" / "raw" / "{group}" / "{sample}.sam"),
    benchmark:
        BENCHMARK_FP / "bowtie2_{sample}_{group}.tsv"
    log:
        LOG_FP / "bowtie2_{sample}_{group}.log",
    params:
        basename=str(DEMIC_FP / "bowtie2" / "{group}" / "contigs"),
    threads: 4
    conda:
        "sbx_demic_bio_env.yml"
    resources:
        mem_mb=8000,
        runtime=240,
    shell:
        """
        bowtie2 -q -x {params.basename} \
        -1 {input.rp1} -2 {input.rp2} -p {threads} \
        -S {output} \
        2>&1 | tee {log}
        """


rule samtools_sort:
    input:
        str(MAPPING_FP / "demic" / "raw" / "{group}" / "{sample}.sam"),
    output:
        temp_files=temp(
            str(MAPPING_FP / "demic" / "sorted" / "{group}" / "{sample}.bam")
        ),
        sorted_files=str(MAPPING_FP / "demic" / "sorted" / "{group}" / "{sample}.sam"),
    log:
        str(MAPPING_FP / "demic" / "logs" / "samtools_{sample}_{group}.log"),
    benchmark:
        BENCHMARK_FP / "samtools_sort_{sample}_{group}.tsv"
    threads: 4
    conda:
        "sbx_demic_bio_env.yml"
    resources:
        mem_mb=8000,
        runtime=240,
    shell:
        """
        samtools view -@ {threads} -bS {input} | \
        samtools sort -@ {threads} - -o {output.temp_files}
        samtools view -@ {threads} -h {output.temp_files} > {output.sorted_files}
        """


rule metabat2:
    input:
        contigs=COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
        bams=expand(MAPPING_FP / "demic" / "sorted" / "{{group}}" / "{sample}.bam", sample=Samples.keys()),
    output:
        directory(DEMIC_FP / "metabat2" / "{group}")
    benchmark:
        BENCHMARK_FP / "metabat2_{group}.tsv"
    log:
        LOG_FP / "metabat2_{group}.log",
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "sbx_demic_bio_env.yml"
    shell:
        """
        mkdir -p {output}
        cd {output}
        runMetaBat.sh {input.contigs} {input.bams}
        """


rule install_demic:
    output:
        out=DEMIC_FP / ".installed",
    conda:
        "sbx_demic_env.yml"
    script:
        "scripts/install_demic.R"


def samples_for_group(group: str) -> list:
    cg = coassembly_groups(Cfg["coassembly_demic"]["group_file"], Samples.keys())
    return [s for g, s in zip(cg[0], cg[1]) if g == group]


def clusters_for_group(group: str) -> list:
    fns = os.listdir(DEMIC_FP / "maxbin")
    return [
        s.split(".")[1]
        for s in fns
        if len(s.split(".")) > 2
        and s.split(".")[0] == group
        and s.split(".")[2] == "fasta"
    ]


def demic_input_for_group(wildcards) -> list:
    checkpoint_output = checkpoints.maxbin.get(**wildcards).output[0]
    return [
        MAPPING_FP
        / "demic"
        / "sorted"
        / f"{wildcards.group}"
        / f"{sample}_{cluster}.sam"
        for sample in samples_for_group(wildcards.group)
        for cluster in clusters_for_group(wildcards.group)
    ]


rule run_demic:
    input:
        #demic_input_for_group,
        expand(
            MAPPING_FP / "demic" / "sorted" / "{{group}}" / "{sample}.sam",
            sample=Samples.keys(),
        ),
        #DEMIC_FP / "maxbin" / "{group}",
        DEMIC_FP / "metabat2" / "{group}",
        DEMIC_FP / ".installed",
    output:
        str(MAPPING_FP / "demic" / "DEMIC_OUT" / "{group}" / "all_PTR.txt"),
    benchmark:
        BENCHMARK_FP / "run_demic_{group}.tsv"
    log:
        LOG_FP / "run_demic_{group}.log",
    params:
        r_installer=get_demic_path() / "envs" / "install.R",
        demic=get_demic_path() / "vendor_demic_v1.0.2" / "DEMIC.pl",
        sam_dir=str(MAPPING_FP / "demic" / "sorted" / "{group}"),
        fasta_dir=str(DEMIC_FP / "maxbin"),
        keep_all=Cfg["sbx_demic"]["keepall"],
        extras=Cfg["sbx_demic"]["extras"],
    threads: 4
    conda:
        "sbx_demic_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        perl {params.demic} --output_all {params.keep_all} {params.extras} \
        --thread_num {threads} \
        -S {params.sam_dir} -F {params.fasta_dir} \
        -O $(dirname {output}) 2>&1 {log}
        """


rule aggregate_demic:
    input:
        expand(
            MAPPING_FP / "demic" / "DEMIC_OUT" / "{group}" / "all_PTR.txt",
            group=list(
                set(
                    coassembly_groups(
                        Cfg["coassembly_demic"]["group_file"], Samples.keys()
                    )[0]
                )
            ),
        ),
    output:
        MAPPING_FP / "demic" / "DEMIC_OUT" / "all_PTR.txt",
    shell:
        "cat {input} > {output}"
