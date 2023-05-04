process REPP_LIST_DB {
    container { params.repp_container }
    label 'process_low'

    input:

    output:
    stdout

    script:
    def repp_repository_env = params.repp_repository
        ? "REPP_DATA_DIR=${params.repp_repository}"
            : ''
    """
    ${repp_repository_env} \
        /go/bin/repp list database
    """
}
