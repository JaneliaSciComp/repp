package main

import (
	"log"
	"os/exec"

	"github.com/Lattice-Automation/repp/internal/cmd"
)

func main() {
	checkDependencies()

	if err := cmd.RootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}

func checkDependencies() {
	if _, err := exec.LookPath("blastn"); err != nil {
		log.Fatal(`No blastn found. Is BLAST installed? https://blast.ncbi.nlm.nih.gov/Blast.cgi`)
	}

	if _, err := exec.LookPath("blastdbcmd"); err != nil {
		log.Fatal(`No blastdbcmd found. Is BLAST installed? https://blast.ncbi.nlm.nih.gov/Blast.cgi`)
	}

	if _, err := exec.LookPath("primer3_core"); err != nil {
		log.Fatal(`No primer3_core found. Is Primer3 installed? https://primer3.org/manual.html`)
	}

	if _, err := exec.LookPath("ntthal"); err != nil {
		log.Fatal(`No ntthal found. Is Primer3 installed? https://primer3.org/manual.html`)
	}
}
