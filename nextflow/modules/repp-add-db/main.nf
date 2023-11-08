process REPP_ADD_DB {
    container { params.repp_container }
    label 'process_low'

    input:
    tuple val(db_name), path(db_path), val(db_cost)

    output:
    tuple val(db_name), path(db_path)

    script:
    def mkdir_repp_repo = params.repp_repository
        ? "mkdir -p ${params.repp_repository}"
        : ''
    def repp_repository_env = params.repp_repository
        ? "REPP_DATA_DIR=${params.repp_repository}"
            : ''
    """
    ${mkdir_repp_repo}
    ${repp_repository_env} \
        /go/bin/repp add database --name ${db_name} -c ${db_cost} ${db_path}
    """
}
