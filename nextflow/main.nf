#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    REPP_ADD_DB;
} from './modules/repp-add-db/main'

include {
    REPP_MAKE_PLASMID;
} from './modules/repp-make-plasmid/main'

workflow {
    main:
    if (params.cmd == "add-db") {
        log.info "!!!! REPP REPO: ${params.repp_repository}"
        REPP_ADD_DB(
            Channel.of(
                [
                    params.dbname,
                    params.dbpath,
                    params.dbcost,
                ]
            )
        )
    } else if (params.cmd == "make-plasmid") {
        REPP_MAKE_PLASMID(
            Channel.of(
                [
                    params.seq_file,
                    params.assembly_output,
                    params.dbs,
                ]
            )
        )
    } else {
        log.error "Invalid command: ${params.cmd}"
    }
}
