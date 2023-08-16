# Softwares used in the biomarker module

## [cnvkit2scarhrd](https://github.com/hydra-genetics/biomarker/blob/develop/workflow/scripts/cnvkit2scarhrd.py)
Script that modifies the cnvkit segmentation file so that it can be used by scarhrd.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__scarhrd__cnvkit2scarhrd#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__scarhrd__cnvkit2scarhrd#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit2scarhrd#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit2scarhrd#

---

## [fix_scarhrd_output](https://github.com/hydra-genetics/biomarker/blob/develop/workflow/scripts/fix_scarhrd_output.py)
Script that modifies the scarhrd output for improved readability of the resulting score.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__scarhrd__fix_scarhrd_output#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__scarhrd__fix_scarhrd_output#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__fix_scarhrd_output#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__fix_scarhrd_output#

---

## [msisensor_pro](https://github.com/xjtu-omics/msisensor-pro)
Calculates the % of microsatelites detected that are instable which can be used to determine MSS or MSI status.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__msisensor_pro__msisensor_pro#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__msisensor_pro__msisensor_pro#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__msisensor_pro#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__msisensor_pro#

---

## [scarhrd](https://github.com/sztup/scarHRD)
Calculates a HRD score by the sum of three different scores that counts the number of cnv events found. The program uses an external segmentation file from for example cnvkit.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__scarhrd__scarhrd#

#### :left_right_arrow: input / output files

### SNAKEMAKE_RULE_TABLE__scarhrd__scarhrd#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__scarhrd#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__scarhrd#

---

## [tmb](https://github.com/hydra-genetics/biomarker/blob/develop/workflow/scripts/tmb.py)
Python script that calculates the tumor mutational burden. Filters the vcf file based on vep annotations to obtain a confident set of somatic nsSNVs. The file can optionally be filtered by an artifact and background noise file. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__tmb__tmb#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__tmb__tmb#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__tmb#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__tmb#

---