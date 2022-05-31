package config

import _ "embed"

//go:embed enzymes.tsv
var enzymes []byte

//go:embed features.tsv
var features []byte