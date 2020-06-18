import glob


rule cellranger_count:
    input:
        fastqs=lambda wildcards: samples["fastqs"][wildcards.sample]
    output:
        directory("results/cellranger_count/{sample}")
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
    log: stderr="results/logs/cellranger_count/{sample}.log"
    # TODO: make --expect-cells optional (builtin default is 3000)
    # TODO: allow --jobmode=sge (in which case --localcores and --localmem do not apply)
    shell:
        """
        cellranger count --id={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={input.fastqs} \
        {params.sample_option} \
        --expect-cells={params.expect_cells} \
        {params.runtime_options} \
        2> {log.stderr} &&
        mv {wildcards.sample} {output}
        """
