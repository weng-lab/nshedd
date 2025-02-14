import pandas as pd
import glob
import os

#samples = os.listdir('/pi/zhiping.weng-umw/data/sheddn/igSCREEN/ChromBPNet/data')
samples = [s.replace('/pi/zhiping.weng-umw/data/sheddn/igSCREEN/ChromBPNet/data/','') for s in glob.glob('/pi/zhiping.weng-umw/data/sheddn/igSCREEN/ChromBPNet/data/*-SRX*') if not s.startswith('/pi/zhiping.weng-umw/data/sheddn/igSCREEN/ChromBPNet/data/SU')]
folds = ["fold_0", "fold_1", "fold_2", "fold_3", "fold_4"]
dir = "/pi/zhiping.weng-umw/data/sheddn/igSCREEN/ChromBPNet"

rule all:
    input:
        expand("data/{sample}/{sample}_peaks_noblacklist.narrowPeak", sample=samples),
        expand("data/{sample}/{sample}_{fold}_negatives.bed", sample=samples, fold=folds),
        expand("results/{sample}/bias/{fold}/models/{sample}_bias.h5", sample=samples, fold=folds),
        expand("results/{sample}/chrombpnet_model/{fold}/models/chrombpnet_nobias.h5", sample=samples, fold=folds),
        expand("results/{sample}/predictions/{fold}/{sample}_{fold}_chrombpnet_nobias.bw", sample=samples, fold=folds),
        expand("results/{sample}/contributions/{fold}/{sample}_{fold}.profile_scores.bw", sample=samples, fold=folds),
        expand("results/{sample}/contributions/{fold}/{sample}_{fold}.counts_scores.bw", sample=samples, fold=folds),
        expand("results/{sample}/contributions/folds_merged/{sample}.profile_scores.bw", sample=samples),
        expand("results/{sample}/contributions/folds_merged/{sample}.counts_scores.bw", sample=samples),
        expand("results/{sample}/predictions/folds_merged/{sample}_chrombpnet_nobias.bw", sample=samples)

rule blacklist:
    input: 'data/{sample}/{sample}_peaks.narrowPeak'
    output: 'data/{sample}/{sample}_peaks_noblacklist.narrowPeak'
    threads: 1
    log: "log/blacklist/{sample}.log"
    resources:
        nodes=1,
        queue='short',
        runtime=120,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process"
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec --nv /home/nicole.shedd-umw/bin/bioinformatics.sif bedtools intersect -v -a {input} -b resources/temp.bed  > {output}
        ) &> {log}
        """

rule prep_nonpeaks:
    input: 
        peaks = 'data/{sample}/{sample}_peaks_noblacklist.narrowPeak',
        fold = 'splits/{fold}.json'
    output: 'data/{sample}/{sample}_{fold}_negatives.bed'
    log: 'log/prep_nonpeaks/{sample}_{fold}.log'
    resources:
        nodes=10,
        queue='gpu',
        runtime=40,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process",
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            trap "rm -rf data/{wildcards.sample}/{wildcards.sample}_{wildcards.fold}_auxiliary" EXIT
            singularity exec --nv ~/bin/chrombpnet_nicole.sif \
                chrombpnet prep nonpeaks \
                -g resources/hg38.fa \
                -p {input.peaks} \
                -c resources/hg38.chrom.sizes \
                -fl {input.fold} \
                -br resources/blacklist.bed.gz \
                -o data/{wildcards.sample}/{wildcards.sample}_{wildcards.fold}
                -st $(( 1000 / (2 ** ({snakemake.attempt} - 1)) > 1 ? 1000 / (2 ** ({snakemake.attempt} - 1)) : 1 ))
        
        ) &> {log}
        """

rule train_bias: ## try to add reading .html file an increasing -b parameter if it errors out and write a file with the final bias parameter value
    input: 
        peaks = "data/{sample}/{sample}_peaks_noblacklist.narrowPeak",
        bam = "data/{sample}/{sample}.bam",
        negatives = "data/{sample}/{sample}_{fold}_negatives.bed",
        fold = "splits/{fold}.json"
    output: "results/{sample}/bias/{fold}/models/{sample}_bias.h5"
    log: 'log/train_bias/{sample}_{fold}.log'
    resources:
        nodes=10,
        queue='gpu',
        runtime=8640,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process"
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            tmpDir=/tmp/sheddn/{wildcards.sample}_{wildcards.fold}/
            mkdir -p $tmpDir; trap "rm -rf $tmpDir" EXIT
            singularity exec --nv ~/bin/chrombpnet_nicole.sif \
                chrombpnet bias pipeline \
                -ibam {input.bam} \
                -d 'ATAC' \
                -g resources/hg38.fa \
                -c resources/hg38.chrom.sizes \
                -p {input.peaks} \
                -n {input.negatives} \
                -fl {input.fold} \
                -b $((0.5 + (0.1 * ({snakemake.attempt} - 1) ))) \
                -o $tmpDir \
                -fp {wildcards.sample}
        mkdir -p results/{wildcards.sample}/bias/{wildcards.fold}/
        cp -r $tmpDir/* results/{wildcards.sample}/bias/{wildcards.fold} 
            
        ) &> {log}
        """

rule train_chrombpnet: ## check quality report in chrombpnet model???
    input: 
        peaks = 'data/{sample}/{sample}_peaks_noblacklist.narrowPeak',
        bam = 'data/{sample}/{sample}.bam',
        negatives = 'data/{sample}/{sample}_{fold}_negatives.bed',
        fold = 'splits/{fold}.json',
        bias = 'results/{sample}/bias/{fold}/models/{sample}_bias.h5'
    output: 'results/{sample}/chrombpnet_model/{fold}/models/chrombpnet_nobias.h5'
    log: 'log/train_chrombpnet/{sample}_{fold}.log'
    resources:
       nodes=10,
        queue='gpu',
        runtime=8640,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process"
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            tmpDir=/tmp/sheddn/{wildcards.sample}_{wildcards.fold}/
            mkdir -p $tmpDir; trap "rm -rf $tmpDir" EXIT
            singularity exec --nv ~/bin/chrombpnet_nicole.sif \
                chrombpnet pipeline \
                -ibam {input.bam} \
                -d 'ATAC' \
                -g resources/hg38.fa \
                -c resources/hg38.chrom.sizes \
                -p {input.peaks} \
                -n {input.negatives} \
                -fl {input.fold} \
                -b {input.bias} \
                -o $tmpDir
            mkdir -p results/{wildcards.sample}/chrombpnet_model/{wildcards.fold}/
            cp -r $tmpDir/* results/{wildcards.sample}/chrombpnet_model/{wildcards.fold}
            
        ) &> {log}
        """


rule prediction_bigwig:
    input: 'results/{sample}/chrombpnet_model/{fold}/models/chrombpnet_nobias.h5'
    output: 'results/{sample}/predictions/{fold}/{sample}_{fold}_chrombpnet_nobias.bw'
    log: 'log/prediction/{sample}_{fold}.log'
    resources:
        nodes=10,
        queue='gpu',
        runtime=8640,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process"
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec --nv ~/bin/chrombpnet_nicole.sif \
                chrombpnet pred_bw \
                -cmb {input} \
                -r resources/iCREs_and_bCREs.narrowPeak \
                -g resources/hg38.fa \
                -c resources/hg38.chrom.subset.sizes \
                -op results/{wildcards.sample}/predictions/{wildcards.fold}/{wildcards.sample}_{wildcards.fold}
        ) &> {log}
        """

rule profile_contribution_bigwig:
    input: 'results/{sample}/chrombpnet_model/{fold}/models/chrombpnet_nobias.h5'
    output: 'results/{sample}/contributions/{fold}/{sample}_{fold}.profile_scores.bw'
    log: 'log/profile_contribution/{sample}_{fold}.log'
    resources:
        nodes=10,
        queue='gpu',
        runtime=8640,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process"
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec --nv ~/bin/chrombpnet_nicole.sif \
                chrombpnet contribs_bw \
                -m {input} \
                -r resources/iCREs_and_bCREs.narrowPeak \
                -g resources/hg38.fa \
                -c resources/hg38.chrom.subset.sizes \
                -op results/{wildcards.sample}/contributions/{wildcards.fold}/{wildcards.sample}_{wildcards.fold} \
                -pc profile
        ) &> {log}
        """

rule counts_contribution_bigwig:
    input: 'results/{sample}/chrombpnet_model/{fold}/models/chrombpnet_nobias.h5'
    output: 'results/{sample}/contributions/{fold}/{sample}_{fold}.counts_scores.bw'
    log: 'log/counts_contribution/{sample}_{fold}.log'
    resources:
        nodes=10,
        queue='gpu',
        runtime=8640,
        mem_mb=10000,
        gpu="num=1:gmodel=TeslaV100_SXM2_32GB:mode=exclusive_process"
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec --nv ~/bin/chrombpnet_nicole.sif \
                chrombpnet contribs_bw \
                -m {input} \
                -r resources/iCREs_and_bCREs.narrowPeak \
                -g resources/hg38.fa \
                -c resources/hg38.chrom.subset.sizes \
                -op results/{wildcards.sample}/contributions/{wildcards.fold}/{wildcards.sample}_{wildcards.fold} \
                -pc counts
        ) &> {log}
        """

rule merge_profile:
    input: 
        fold_0='results/{sample}/contributions/fold_0/{sample}_fold_0.profile_scores.bw',
        fold_1='results/{sample}/contributions/fold_1/{sample}_fold_1.profile_scores.bw',
        fold_2='results/{sample}/contributions/fold_2/{sample}_fold_2.profile_scores.bw',
        fold_3='results/{sample}/contributions/fold_3/{sample}_fold_3.profile_scores.bw',
        fold_4='results/{sample}/contributions/fold_4/{sample}_fold_4.profile_scores.bw',
    output: 'results/{sample}/contributions/folds_merged/{sample}.profile_scores.bw'
    log: 'log/merge_profile/{sample}.log'
    resources:
        nodes=1,
        queue='short',
        runtime=480,
        mem_mb=10000
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec /home/nicole.shedd-umw/bin/bioinformatics.sif bigwigAverage -b {input.fold_0} {input.fold_1} {input.fold_2} {input.fold_3} {input.fold_4} -o {output}
        ) &> {log}
        """

rule merge_counts:
    input: 
        fold_0='results/{sample}/contributions/fold_0/{sample}_fold_0.counts_scores.bw',
        fold_1='results/{sample}/contributions/fold_1/{sample}_fold_1.counts_scores.bw',
        fold_2='results/{sample}/contributions/fold_2/{sample}_fold_2.counts_scores.bw',
        fold_3='results/{sample}/contributions/fold_3/{sample}_fold_3.counts_scores.bw',
        fold_4='results/{sample}/contributions/fold_4/{sample}_fold_4.counts_scores.bw',
    output: 'results/{sample}/contributions/folds_merged/{sample}.counts_scores.bw'
    log: 'log/merge_counts/{sample}.log'
    resources:
        nodes=1,
        queue='short',
        runtime=480,
        mem_mb=10000
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec /home/nicole.shedd-umw/bin/bioinformatics.sif bigwigAverage -b {input.fold_0} {input.fold_1} {input.fold_2} {input.fold_3} {input.fold_4} -o {output}
        ) &> {log}
        """

rule merge_predictions:
    input: 
        fold_0='results/{sample}/predictions/fold_0/{sample}_fold_0_chrombpnet_nobias.bw',
        fold_1='results/{sample}/predictions/fold_1/{sample}_fold_1_chrombpnet_nobias.bw',
        fold_2='results/{sample}/predictions/fold_2/{sample}_fold_2_chrombpnet_nobias.bw',
        fold_3='results/{sample}/predictions/fold_3/{sample}_fold_3_chrombpnet_nobias.bw',
        fold_4='results/{sample}/predictions/fold_4/{sample}_fold_4_chrombpnet_nobias.bw',
    output: 'results/{sample}/predictions/folds_merged/{sample}_chrombpnet_nobias.bw'
    log: 'log/merge_predictions/{sample}.log'
    resources:
        nodes=1,
        queue='short',
        runtime=480,
        mem_mb=10000
    params:
        gpu_flag=lambda wildcards, resources: f"-gpu {resources.gpu}" if resources.queue == "gpu" else ""
    shell:
        """
        (
            singularity exec /home/nicole.shedd-umw/bin/bioinformatics.sif bigwigAverage -b {input.fold_0} {input.fold_1} {input.fold_2} {input.fold_3} {input.fold_4} -o {output}
        ) &> {log}
        """

