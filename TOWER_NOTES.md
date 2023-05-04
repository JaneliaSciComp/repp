### Setup your tower environment

Follow the [instructions from the wiki](https://wikis.janelia.org/display/SCSW/Using+Nextflow+Tower) to setup your tower environment.

For the pipeline to launch use `https://github.com/JaneliaSciComp/repp` and select `lsf` in your *Config Profiles* option.

### Launch repp job

Once your tower environment is setup launch the pipeline from the launchpad. You have to select the command and fill in the parameters from the corresponding command group (`Add database options` or `Make plasmid pptions`). The examples folder also contains two JSON files one for [adding a database](examples/add-db.json) and one for [running the assembly](examples/make-plasmid.json). If you use the default CSV output format, and you don't specify anything for the output, the result CSV files will be generated in the same folder as the input and their names will be `<inputfilenamewithoutextension>-strategy.csv` and `<inputfilenamewithoutextension>-reagents.csv`
