process REPP_MAKE_PLASMID {
    container { params.repp_container }
    cpus { params.make_plasmid_cpus }

    input:
    tuple path(input_seq),
          path(oligos_manifest),
	  path(output_path),
	  val(db_names)

    output:
    tuple path(input_seq)

    script:
    def repp_repository_env = params.repp_repository
    	? "REPP_DATA_DIR ${params.repp_repository}"
	: ''
    """
    ${repp_repository_env} \
        /app/bin/repp make sequence \
	-i ${input_seq} \
	${dbs_arg} \
	${oligos_manifest_arg} \
	${output_arg} \
	${output_format_arg}

    """
}
