process REPP_ADD_DB {
    container { params.repp_container }
    label 'process_low'

    input:
    tuple val(db_name), path(db_path), val(db_cost)

    output:
    tuple val(db_name), path(db_path)

    script:
    def repp_repository_env = params.repp_repository
    	? "REPP_DATA_DIR ${params.repp_repository}"
	: ''
    """
    ${repp_repository_env} \
        /app/bin/repp add database --name ${db_name} -c ${db_cost} ${db_path}
    """
}
