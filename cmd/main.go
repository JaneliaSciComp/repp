package main

import (
	"log"

	"github.com/Lattice-Automation/repp/internal/cmd"
)

func main() {
	if err := cmd.RootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}
