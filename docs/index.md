# Hydra-genetics

We are an organization/community with the goal of making [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline development easier, faster, a bit more structured and of higher quality.

We do this by providing [snakemake modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) that can be combined to create a complete analysis or included in already existing pipelines. All modules are subjected to extensive testing to make sure that new releases doesn't unexpectedly break existing pipeline or deviate from guidelines and best practices on how to write code.

# biomarker module
The biomarker module consists of programs used for producing biomarker values. The biomarkers currently implemented are:

 - TMB: tumor mutational burden
 - MSI: micro satellite instability score
 - **Under development** HRD: homologous recombination deﬁciency score
