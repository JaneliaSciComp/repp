include {
    get_runtime_opts;
    normalized_file_name;
    parentfile;
} from '../../lib/utils'

process REPP_MAKE_PLASMID {
    container { params.repp_container }
    containerOptions { get_runtime_opts([
        parentfile(input_seq, 1),
        parentfile(params.primers_databases, 1),
        parentfile(params.synth_frags_databases, 1)
        parentfile(assembly_output, 2),
    ]) }
    cpus { params.make_plasmid_cpus }
    memory { "${params.make_plasmid_mem_gb} GB"}

    input:
    tuple val(input_seq),
          val(assembly_output),
          val(db_names)

    output:
    val(input_seq)

    script:
    def repp_repository_env = params.repp_repository
        ? "REPP_DATA_DIR=${params.repp_repository}"
            : ''
    def dbs_arg = db_names
        ? "--dbs ${db_names}"
        : ''
    def primers_databases_arg = params.primers_databases
        ? "-m ${normalized_file_name(params.primers_databases)}"
        : ''
    def synth_frags_databases_arg = params.synth_frags_databases
        ? "-s ${params.synth_frags_databases}"
        : ''
    def max_solutions_arg = params.max_solutions
        ? "-n ${params.max_solutions}"
        : ''
    def config_arg = params.config
        ? "--config ${params.config}"
        : ''
    def assembly_output_name = assembly_output
    def mk_output_dir
    def output_arg
    if (assembly_output_name) {
        def assembly_output_dir = file(assembly_output_name).parent
        mk_output_dir = "mkdir -p ${assembly_output_dir}"
        output_arg = "-o ${normalized_file_name(assembly_output_name)}"
    } else {
        def input_dir = file(input_seq).parent
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
            -i ${normalized_file_name(input_seq)} \
            ${dbs_arg} \
            ${primers_databases_arg} \
            ${synth_frags_databases_arg} \
            ${output_arg} \
            ${sequence_identity_arg} \
            ${max_solutions_arg} \
            ${config_arg} \
            ${output_format_arg}
    """
}
