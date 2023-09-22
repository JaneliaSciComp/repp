package main

import (
	"log"
	"os"
	"os/exec"

	"github.com/Lattice-Automation/repp/internal/cmd"
)

func main() {
	checkDependencies()
	if err := cmd.RootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

func checkDependencies() {
	if _, err := exec.LookPath(getExecutable("NCBITOOLS_HOME", "bin", "blastn")); err != nil {
		log.Fatal(`No blastn found. Is BLAST installed? https://blast.ncbi.nlm.nih.gov/Blast.cgi`)
	}

	if _, err := exec.LookPath(getExecutable("NCBITOOLS_HOME", "bin", "blastdbcmd")); err != nil {
		log.Fatal(`No blastdbcmd found. Is BLAST installed? https://blast.ncbi.nlm.nih.gov/Blast.cgi`)
	}

	if _, err := exec.LookPath(getExecutable("NCBITOOLS_HOME", "bin", "makeblastdb")); err != nil {
		log.Fatal(`No makeblastdb found. Is BLAST installed? https://blast.ncbi.nlm.nih.gov/Blast.cgi`)
	}

	if _, err := exec.LookPath(getExecutable("PRIMER3_HOME", "bin", "primer3_core")); err != nil {
		log.Fatal(`No primer3_core found. Is Primer3 installed? https://primer3.org/manual.html`)
	}

	if _, err := exec.LookPath(getExecutable("PRIMER3_HOME", "bin", "ntthal")); err != nil {
		log.Fatal(`No ntthal found. Is Primer3 installed? https://primer3.org/manual.html`)
	}
}

func getExecutable(exeHomeEnvVar, binSubDir, exeName string) string {
	exeHome := os.Getenv(exeHomeEnvVar)
	if exeHome == "" {
		// if no home or install dir is set, assume it's in the PATH
		return exeName
	}
	if binSubDir == "" {
		return exeHome + "/" + exeName
	} else {
		return exeHome + "/" + binSubDir + "/" + exeName
	}
}
