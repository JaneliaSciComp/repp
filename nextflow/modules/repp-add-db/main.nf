process REPP_ADD_DB {
    container { params.repp_container }
    label 'process_low'

    input:
    tuple val(db_name), path(db_path), val(db_cost)
    tuple path(repp_repo_parent), val(repp_repo_name)

    output:
    tuple val(db_name), path(db_path)

    script:
    def repp_repository="${repp_repo_parent}/${repp_repo_name}"
    def repp_repository_env = "REPP_DATA_DIR=${repp_repository}"
    """
    repp_repo_parent_fullpath=\$(readlink ${repp_repo_parent})
    repp_repo_fullpath="\${repp_repo_parent_fullpath}/${repp_repo_name}"
    echo "Repp repo dir: \${repp_repo_fullpath}"
    mkdir -p \${repp_repo_fullpath}
    ${repp_repository_env} \
        /go/bin/repp add database --name ${db_name} -c ${db_cost} ${db_path}
    """
}
