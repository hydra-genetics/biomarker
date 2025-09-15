__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("7.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


def compile_output_list(wildcards):
    platform = units.platform.iloc[0]
    output_files = []
    files = {
        "biomarker/msisensor_pro": [""],
        "biomarker/tmb": [".TMB.txt"],
        "biomarker/finaletoolkit_end_motifs": [".end-motifs.tsv"],
        "biomarker/finaletoolkit_interval_end_motifs": [".interval-end-motifs.tsv"],
        "biomarker/finaletoolkit_mds": [".mds.txt"],
        "biomarker/finaletoolkit_interval_mds": [".interval-mds.txt"],
        "biomarker/finaletoolkit_frag_length_bins": [".frag-length-bins.tsv"],
        "biomarker/finaletoolkit_interval_mds": [".interval-mds.txt"],
        "biomarker/fragmentomics_fragment_length_patient_score": [".fragment_length_patient_score.txt"],
    }
    output_files += [
        f"{prefix}/{sample}_{unit_type}{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for platform in units.loc[(sample,)].platform
        if platform not in ["ONT", "PACBIO"]
        for unit_type in get_unit_types(units, sample)
        if unit_type in ["N", "T"]
        for suffix in files[prefix]
    ]
    return output_files
