process REPP_ADD_DB {
    container { params.repp_container }
    label 'process_low'

    input:
    tuple val(db_name), path(db_path), val(db_cost)
    path(repp_repository)

    output:
    tuple val(db_name), path(db_path)

    script:
    def repp_repository_env = repp_repository
        ? "REPP_DATA_DIR=${repp_repository}"
            : ''
    """
    repp_repo_fullpath=\$(realpath ${repp_repository})
    echo "Repp repo dir: \${repp_repo_fullpath}"
    mkdir -p \${repp_repo_fullpath}
    ${repp_repository_env} \
        /go/bin/repp add database --name ${db_name} -c ${db_cost} ${db_path}
    """
}
