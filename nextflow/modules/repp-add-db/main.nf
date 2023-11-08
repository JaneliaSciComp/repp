include {
    get_runtime_opts;
    parentfile;
} from '../../lib/utils'

process REPP_ADD_DB {
    container { params.repp_container }
    containerOptions { get_runtime_opts([
        parentfile(repp_repo_parent, 1),
    ]) }
    label 'process_low'

    input:
    tuple path(db_path), val(db_name), val(db_cost)
    tuple path(repp_repo_parent), val(repp_repo_name)

    output:
    tuple val(db_name), path(db_path), env(repp_repo_fullpath)

    script:
    def repp_repository="${repp_repo_parent}/${repp_repo_name}"
    def repp_repository_env = "REPP_DATA_DIR=${repp_repository}"
    """
    umask 0002
    repp_repo_parent_fullpath=\$(readlink ${repp_repo_parent})
    repp_repo_fullpath="\${repp_repo_parent_fullpath}/${repp_repo_name}"
    echo "Repp repo dir: \${repp_repo_fullpath}"
    mkdir -p \${repp_repo_fullpath}
    ${repp_repository_env} \
        /go/bin/repp add database --name ${db_name} -c ${db_cost} ${db_path}
    """
}
