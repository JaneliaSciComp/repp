include {
    get_runtime_opts;
    normalized_file_name;
    parentfile;
} from '../../lib/utils'

process REPP_MAKE_PLASMID {
    container { params.repp_container }
    containerOptions { get_runtime_opts([
        parentfile(params.oligos_manifest, 1),
        parentfile(assembly_output, 2),
    ]) }
    cpus { params.make_plasmid_cpus }
    memory { "${params.make_plasmid_mem_gb} GB"}

    input:
    tuple path(input_seq),
	      val(assembly_output),
	      val(db_names)

    output:
    path(input_seq)

    script:
    def repp_repository_env = params.repp_repository
    	? "REPP_DATA_DIR=${params.repp_repository}"
	    : ''
    def dbs_arg = db_names
        ? "--dbs ${db_names}"
        : ''
    def oligos_manifest_arg = params.oligos_manifest
        ? "-m ${normalized_file_name(params.oligos_manifest)}"
        : ''
    def mk_output_dir
    def output_arg
    if (assembly_output) {
        assembly_output_dir = file(assembly_output).parent
        mk_output_dir = "mkdir -p ${assembly_output_dir}"
        output_arg = "-o ${normalized_file_name(assembly_output)}"
    } else {
        mk_output_dir = ''
        output_arg = ''
    }
    def output_format_arg = params.assembly_output_format
        ? "-f ${params.assembly_output_format}"
        : ''
    def sequence_identity_arg = params.sequence_identity > 0 && params.sequence_identity < 100
        ? "-p ${params.sequence_identity}"
        : ''
    """
    ${mk_output_dir}

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
