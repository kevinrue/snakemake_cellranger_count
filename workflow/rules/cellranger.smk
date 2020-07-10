from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

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
        local_memory = get_local_memory(),
        option_str += f"--jobmode=local --localcores={threads} --localmem={local_memory}"
    elif jobmode == "sge":
        memory_per_cpu = config['cellranger']['memory_per_cpu']
        option_str += f"--jobmode=sge --maxjobs={threads} --mempercore={memory_per_cpu}"
    else:
        raise NameError(f"Invalid job mode: {jobmode}")
    
    return option_str


def get_local_memory():
    '''
    Get cellranger mkref memory options.
    '''
    jobmode = config['cellranger']['jobmode']
    threads = config['cellranger']['threads']
    
    memory_value = config['cellranger']['memory_per_cpu'] * threads
    
    return memory_value



def get_rule_threads():
    '''
    Get the number of threads given to run the rule.
    '''
    threads = 1
    
    jobmode = config['cellranger']['jobmode']
    if jobmode == "local":
        threads = config['cellranger']['threads']
    elif jobmode == "sge":
        threads = 1
    else:
        raise NameError(f"Invalid job mode: {jobmode}")
    
    return threads


rule genome:
    input:
        FTP.remote(config['cellranger']['genome'])
    output:
        'resources/genome.fa.gz'
    shell:
        '''
        mv {input} {output}
        '''


rule genesets:
    input:
        FTP.remote(config['cellranger']['genesets'])
    output:
        'resources/genesets.gtf'
    shell:
        '''
        mv {input} {output}
        '''
        

rule cellranger_mkref:
    input:
        fasta='resources/genome.fa.gz',
        genes='resources/genesets.gtf'
    output:
        directory("resources/cellranger_index")
    params:
        memory=get_local_memory(),
        threads=config['cellranger']['threads']
    envmodules:
        "bio/cellranger/3.1.0"
    threads: config['cellranger']['threads']
    resources:
        mem_free_gb=config['cellranger']['memory_per_cpu']
    log:
        err="results/logs/cellranger_mkref.err",
        out="results/logs/cellranger_mkref.out",
    shell:
        """
        gunzip -c {input.genes} > {input.genes}.tmp &&
        gunzip -c {input.fasta} > {input.fasta}.tmp &&
        cellranger mkref \
        --genome=cellranger_index \
        --fasta="{input.fasta}.tmp" \
        --genes="{input.genes}.tmp" \
        --nthreads={params.threads} \
        --memgb={params.memory}
        2> {log.err} > {log.out} &&
        mv cellranger_index {output} &&
        rm {input.genes}.tmp {input.fasta}.tmp
        """


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
        "bio/cellranger/3.1.0"
    threads: get_rule_threads()
    resources:
        mem_free_gb=config['cellranger']['memory_per_cpu']
    log:
        err="results/logs/cellranger_count/{sample}.err",
        out="results/logs/cellranger_count/{sample}.out"
    shell:
        """
        cellranger count --id={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={input.fastqs} \
        {params.sample_option} \
        --expect-cells={params.expect_cells} \
        {params.runtime_options} \
        2> {log.err} > {log.out} &&
        rm -rf results/cellranger_count/{wildcards.sample} &&
        mv {wildcards.sample} results/cellranger_count/
        """
