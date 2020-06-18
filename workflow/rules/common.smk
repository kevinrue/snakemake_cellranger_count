from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

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
    print(sample_prefix)

    if sample_prefix != ".":
        option_str += f"--sample={sample_prefix}"
        
    return option_str

