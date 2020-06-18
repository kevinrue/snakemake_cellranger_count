from snakemake.utils import validate
import pandas as pd

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

##### helper functions #####

def get_sample_option(wildcards):
    '''
    Build the optional --sample option for `cellranger count`.
    '''
    
    option_str = ""
    
    sample_prefix = samples['prefix'][wildcards.sample]

    if sample_prefix != ".":
        option_str += f"--sample={sample_prefix}"
        
    return option_str


def get_runtime_options():
    '''
    
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

