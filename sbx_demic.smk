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
            rp=coassembly_groups(Cfg["coassembly_demic"]["group_file"], Samples.keys())[
                2
            ],
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
        LOG_FP / "prep_samples_for_concatenation_paired_{sample}_{group}.log",
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
        rules.all_coassemble.input.b,
    output:
        r1=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz"),
        r2=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz"),
    params:
        w1=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_1.fastq")),
        w2=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_2.fastq")),
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "sbx_coassembly_env.yml"
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
        LOG_FP / "coassemble_paired_{group}.log",
    params:
        assembly_dir=str(COASSEMBLY_DEMIC_FP / "{group}"),
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "sbx_coassembly_env.yml"
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
        DEMIC_FP / "decontam_list.txt",
    params:
        decontam_fp=QC_FP / "decontam",
    shell:
        "find {params.decontam_fp} -iname '*.fastq.gz' > {output}"


rule maxbin_download:
    """Supposedly bioconda-recipes has updated the maxbin recipe to not point
    at the JBEI source (which is down) but no version of conda/bioconda I've 
    tried works
    https://github.com/bioconda/bioconda-recipes/pull/25064
    """
    output:
        get_demic_path() / ".maxbin.installed"
    params:
        target_dir=get_demic_path()
    conda:
        "demic_bio_env.yml"
    shell:
        """
        ENV_DIR=$(dirname $(dirname $(which bowtie2)))
        wget -qO- https://anaconda.org/bioconda/maxbin2/2.2.7/download/linux-64/maxbin2-2.2.7-h87f3376_4.tar.bz2 \
        | tar xj -C $ENV_DIR
        touch {output}
        """


rule maxbin_setup:
    """Because everything sucks, the tarball for MaxBin has absolute paths to 
    the owner's computer
    """
    input:
        get_demic_path() / "MaxBin-2.2.7" / "setting"
    output:
        get_demic_path() / "MaxBin-2.2.7" / ".fixed"
    conda:
        "demic_bio_env.yml"
    shell:
        """
        rm {input}
        echo "[FragGeneScan] $(which FragGeneScan)" >> {input}
        echo "[Bowtie2] $(which bowtie2)" >> {input}
        echo "[HMMER3] $(which hmmscan)" >> {input}
        touch {output}
        """


rule maxbin:
    input:
        a=expand(
            COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa",
            group=list(
                set(
                    coassembly_groups(
                        Cfg["coassembly_demic"]["group_file"], Samples.keys()
                    )[0]
                )
            ),
        ),
        b=rules.all_coassemble.input.b,
        decontam_list=DEMIC_FP / "decontam_list.txt",
    output:
        DEMIC_FP / "maxbin" / "all_final_contigs.fa",
    benchmark:
        BENCHMARK_FP / "maxbin.tsv"
    log:
        LOG_FP / "maxbin.log",
    params:
        contigs_fasta=str(COASSEMBLY_DEMIC_FP / "all_final_contigs.fa"),
        out_dir=str(DEMIC_FP / "maxbin"),
    conda:
        "demic_bio_env.yml"
    shell:
        """
        conda install -c bioconda perl-lwp-simple
        run_MaxBin.pl -thread 10 -contig {params.contigs_fasta} \
        -out {params.out_dir} -reads {input.decontam_list} \
        -verbose 2>&1 | tee {log}
        """


rule bowtie2_build:
    input:
        DEMIC_FP / "maxbin" / "all_final_contigs.fa",
    output:
        [
            DEMIC_FP / "bowtie" / ("contigs" + ext)
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
        BENCHMARK_FP / "bowtie2_build.tsv"
    log:
        LOG_FP / "bowtie2-build.log",
    params:
        basename=str(DEMIC_FP / "bowtie2" / "contigs"),
    threads: 4
    conda:
        "demic_bio_env.yml"
    shell:
        "bowtie2-build --threads {threads} {input} {params.basename} 2>&1 | tee {log}"


# Run bowtie2 with index
rule bowtie2:
    input:
        rules.bowtie2_build.output,
        rp1=expand(
            str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
            sample=Samples.keys(),
        ),
        rp2=expand(
            str(QC_FP / "decontam" / "{sample}_2.fastq.gz"),
            sample=Samples.keys(),
        ),
    output:
        temp(MAPPING_FP / "demic" / "raw" / "{sample}.sam"),
    benchmark:
        BENCHMARK_FP / "bowtie2_{sample}.tsv"
    log:
        LOG_FP / "bowtie2_{sample}.log",
    params:
        basename=str(DEMIC_FP / "bowtie2" / "contigs"),
    threads: 4
    conda:
        "demic_bio_env.yml"
    shell:
        """
        bowtie2 -q -x {params.basename} \
        -1 {input.rp1} -2 {input.rp2} -p {threads} \
        -S {output} \
        2>&1 | tee {log}
        """


rule samtools_sort:
    input:
        str(MAPPING_FP / "demic" / "raw" / "{sample}.sam"),
    output:
        temp_files=temp(str(MAPPING_FP / "demic" / "sorted" / "{sample}.bam")),
        sorted_files=str(MAPPING_FP / "demic" / "sorted" / "{sample}.sam"),
    log:
        str(MAPPING_FP / "demic" / "logs" / "samtools_{sample}.error"),
    benchmark:
        BENCHMARK_FP / "samtools_sort_{sample}.tsv"
    threads: 4
    conda:
        "demic_bio_env.yml"
    shell:
        """
        samtools view -@ {threads} -bS {input} | \
        samtools sort -@ {threads} - -o {output.temp_files} 2>&1 | tee {log}
        samtools view -@ {threads} -h {output.temp_files} > {output.sorted_files} 2>> | tee {log}
        """


rule run_demic:
    input:
        expand(
            str(MAPPING_FP / "demic" / "sorted" / "{sample}.sam"),
            sample=Samples.keys(),
        ),
    output:
        str(MAPPING_FP / "demic" / "DEMIC_OUT" / "all_PTR.txt"),
    benchmark:
        BENCHMARK_FP / "run_demic.tsv"
    log:
        LOG_FP / "run_demic.log",
    params:
        r_installer=get_demic_path() / "envs" / "install.R",
        demic=get_demic_path() / "vendor_demic_v1.0.2" / "DEMIC.pl",
        sam_dir=str(MAPPING_FP / "demic" / "sorted"),
        fasta_dir=str(DEMIC_FP / "maxbin"),
        keep_all=Cfg["sbx_demic"]["keepall"],
        extras=Cfg["sbx_demic"]["extras"],
    threads: 4
    conda:
        "demic_env.yml"
    shell:
        """
        {params.demic} --output_all {params.keep_all} {params.extras} \
        --thread_num {threads} \
        -S {params.sam_dir} -F {params.fasta_dir} \
        -O $(dirname {output}) 2>&1 {log}
        """
