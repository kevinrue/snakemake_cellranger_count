import glob

def get_fastqs(wildcards):
    '''
    Identify pairs of FASTQ files from the sample sheet.
    
    wildcards
    - sample: name of the sample to process.
    '''
    fastqs_glob = f"data/{wildcards.sample}_fastqs"
    fastqs = glob.glob(fastqs_glob)
    
    if len(fastqs) == 0:
        raise OSError(f"No directory matched pattern: {fastqs_glob}")
        
    return {'fastqs' : fastqs}


rule cellranger_count:
    input:
        unpack(get_fastqs)
    output:
        directory("results/cellranger_count/{sample}")
    params:
        transcriptome=config['cellranger']['transcriptome'],
        expect_cells=lambda wildcards, input: samples['expect_cells'][wildcards.sample],
        local_memory=config['cellranger']['memory_per_cpu'] * config['cellranger']['threads'],
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
        --expect-cells={params.expect_cells} \
        --jobmode=local \
        --localcores={params.threads} \
        --localmem={params.local_memory} \
        2> {log.stderr} &&
        mv {wildcards.sample} {output}
        """
