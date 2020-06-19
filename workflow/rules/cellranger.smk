import glob


def get_sample_option(wildcards):
    '''
    Get cellranger sample option.
    '''
    
    option_str = ""
    
    sample_prefix = samples['prefix'][wildcards.sample]

    if sample_prefix != ".":
        option_str += f"--sample={sample_prefix}"
        
    return option_str


def get_runtime_options():
    '''
    Get cellranger runtime options.
    '''
    option_str = ""
    jobmode = config['cellranger']['jobmode']
    threads = config['cellranger']['threads']
    if jobmode == "local":
        local_memory = config['cellranger']['memory_per_cpu'] * threads
        option_str += f"--jobmode=local --localcores={threads} --localmem={local_memory}"
    elif jobmode == "sge":
        memory_per_cpu = config['cellranger']['memory_per_cpu']
        option_str += f"--jobmode=sge --maxjobs={threads} --mempercore={memory_per_cpu}"
    else:
        raise NameError(f"Invalid job mode: {params.jobmode}")
    return option_str


rule cellranger_count:
    input:
        fastqs=lambda wildcards: samples["fastqs"][wildcards.sample]
    output:
        raw="results/cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5",
        filtered="results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        bam="results/cellranger_count/{sample}/outs/possorted_genome_bam.bam"
    params:
        transcriptome=config['cellranger']['transcriptome'],
        expect_cells=lambda wildcards, input: samples['expect_cells'][wildcards.sample],
        runtime_options=get_runtime_options(),
        sample_option=get_sample_option,
        threads=config['cellranger']['threads']
    envmodules:
        "bio/cellranger/3.0.1"
    threads: config['cellranger']['threads']
    resources:
        mem_free_gb=config['cellranger']['memory_per_cpu']
    log: "results/logs/cellranger_count/{sample}.err"
    shell:
        """
        cellranger count --id={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={input.fastqs} \
        {params.sample_option} \
        --expect-cells={params.expect_cells} \
        {params.runtime_options} \
        2> {log} &&
        mv {wildcards.sample} results/cellranger_count/{wildcards.sample}
        """
