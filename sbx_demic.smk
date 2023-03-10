# -*- mode: Snakemake -*-

import os
import sys
from pathlib import Path

COASSEMBLY_FP = ASSEMBLY_FP / "coassembly"
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


def get_demic_path() -> str:
    demic_path = Path(sunbeam_dir) / "extensions" / "sbx_demic"
    if demic_path.exists():
        return demic_path
    raise Error(
        "Filepath for demic not found, are you sure it's installed under extensions/sbx_demic?"
    )


rule all_demic:
    input:
        TARGET_DEMIC,


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


rule maxbin:
    input:
        a=expand(
            COASSEMBLY_FP / "{group}_final_contigs.fa",
            group=list(
                set(
                    coassembly_groups(
                        Cfg["sbx_coassembly"]["group_file"], Samples.keys()
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
        contigs_fasta=str(COASSEMBLY_FP / "all_final_contigs.fa"),
        out_dir=str(DEMIC_FP / "maxbin"),
    conda:
        "demic_bio_env.yml"
    shell:
        """
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
