def get_runtime_opts(paths) {
    def bind_paths = get_mounted_vols_opts(paths)
    "${params.runtime_opts} ${bind_paths}"
}

def parentfile(f, uplevels) {
    if (f) {
        def ff = file(f)
        def parent = ff
        def nlevels = uplevels
        while (nlevels-- > 0) {
            def p = parent.parent;
            if (!p)
                return parent
            else
                parent = p
        
        }
        return parent        
    } else {
        null
    }
}

def get_mounted_vols_opts(paths) {
    def unique_paths = paths
                        .collect { normalized_file_name(it) }
                        .findAll { it ? true : false }
                        .unique(false)
    switch (workflow.containerEngine) {
        case 'docker': 
            return unique_paths
                .collect { "-v $it:$it" }
                .join(' ')
        case 'singularity':
            return unique_paths
                .collect { "-B $it" }
                .join(' ')
        default:
            log.error "Unsupported container engine: ${workflow.containerEngine}"
            ''
    }
}

def normalized_file_name(f) {
    if (f) {
        "${file(f)}"
    } else {
        ''
    }
}
