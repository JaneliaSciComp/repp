process REPP_MAKE_PLASMID {
    container { params.repp_container }
    cpus { params.make_plasmid_cpus }
    memory { "${params.make_plasmid_mem_gb} GB"}

    input:
    tuple path(input_seq),
	      path(assembly_output),
	      val(db_names)

    output:
    tuple path(input_seq)

    script:
    def repp_repository_env = params.repp_repository
    	? "REPP_DATA_DIR=${params.repp_repository}"
	    : ''
    def dbs_arg = db_names
        ? "--dbs ${db_names}"
        : ''
    def oligos_manifest_arg = params.oligos_manifest
        ? "-m ${params.oligos_manifest}"
        : ''
    def output_arg = assembly_output
        ? "-o ${assembly_output}"
        : ''
    def output_format_arg = params.assembly_output_format
        ? "-f ${params.plasmid_assembly_output_format}"
        : ''
    def sequence_identity_arg = params.sequence_identity > 0 && params.sequence_identity < 100
        ? "-p ${params.sequence_identity}"
        : ''
    """
    ${repp_repository_env} \
        /go/bin/repp make sequence \
	-i ${input_seq} \
	${dbs_arg} \
	${oligos_manifest_arg} \
	${output_arg} \
    ${sequence_identity_arg} \
	${output_format_arg}
    """
}
