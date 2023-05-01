${BLAST_PATH} = "C:\Program Files\NCBI\blast-2.14.0+"
${PRIMER3_PATH} = "${env:USERPROFILE}\primer3-latest"
$env:Path = "$env:Path;" + "${BLAST_PATH}\bin;" + "${PRIMER3_PATH}\src;"
