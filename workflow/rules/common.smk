# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

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

if os.path.isfile("config/config.yaml"):

    configfile: "config/config.yaml"


elif os.path.isfile("config.yaml"):

    configfile: "config.yaml"


elif not workflow.overwrite_configfiles:
    raise FileExistsError("No config file found in working directory or passed as argument!")


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
    of = ["biomarker/msisensor_pro/%s_%s" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    of.append(
        ["biomarker/tmb/%s_%s.TMB.txt" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    )
    of.append(
        [
            "biomarker/hrd/%s_%s.hrd_score.txt" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    of.append(
        [
            "biomarker/optitype/%s_%s/%s_%s_hla_type_result.tsv" % (sample, t, sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    return of
