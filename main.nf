#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    REPP_ADD_DB;
} from './nextflow/modules/repp-add-db/main'

include {
    REPP_LIST_DB;
} from './nextflow/modules/repp-list-db/main'

include {
    REPP_MAKE_PLASMID;
} from './nextflow/modules/repp-make-plasmid/main'

workflow {
    main:
    def repp_repo = file(params.repp_repository)

    if (params.reppcmd == "add-db") {
        REPP_ADD_DB(
            Channel.of(
                [
                    params.dbname,
                    params.dbpath,
                    params.dbcost,
                ]
            ),
            [ repp_repo.parent, repp_repo.name ],
        )
    } else if (params.reppcmd == "list-db") {
        REPP_LIST_DB | view
    } else if (params.reppcmd == "make-plasmid") {
        REPP_MAKE_PLASMID(
            Channel.of(
                [
                    params.seq_file,
                    params.assembly_output,
                    params.dbs,
                ]
            ),
            params.make_plasmid_cpus,
            params.make_plasmid_mem_gb,
        )
    } else {
        log.error "Invalid command: ${params.reppcmd}"
    }
}
