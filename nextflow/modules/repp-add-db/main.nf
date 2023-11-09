include {
    get_runtime_opts;
    parentfile;
} from '../../lib/utils'

process REPP_ADD_DB {
    container { params.repp_container }
    containerOptions { 
        get_runtime_opts([
        parentfile(repp_repo_dir, 1),
    ]) }
    label 'process_low'

    input:
    tuple path(db_path), val(db_name), val(db_cost)
    tuple val(repp_repo_dir)

    output:
    tuple val(db_name), path(db_path), val(repp_repo_dir)

    script:
    def mk_repp_repo = repp_repo_dir
                        ? "mkdir -p ${repp_repo_dir}"
                        : ''
    def repp_repository_env = repp_repo_dir
                        ? "REPP_DATA_DIR=${repp_repo_dir}"
                        : ''
    """
    umask 0002
    ${mk_repp_repo}
    ${repp_repository_env} /go/bin/repp add database --name ${db_name} -c ${db_cost} ${db_path}
    """
}
