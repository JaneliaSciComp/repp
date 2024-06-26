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

    if (params.reppcmd == "add-db") {
        REPP_ADD_DB(
            Channel.of(
                [
                    params.dbpath,
                    params.dbname,
                    params.dbcost,
                ]
            ),
            params.repp_repository,
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
            params.repp_repository,
            params.primers_databases ? params.primers_databases : [],
            params.synth_frags_databases ? params.synth_frags_databases : [],
            params.make_plasmid_cpus,
            params.make_plasmid_mem_gb,
        )
    } else {
        log.error "Invalid command: ${params.reppcmd}"
    }
}
