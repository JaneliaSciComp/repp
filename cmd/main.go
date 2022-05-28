package main

import (
	"Lattice-Automation/repp/internal/cmd"
	"log"
)

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func main() {
	if err := cmd.RootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}
