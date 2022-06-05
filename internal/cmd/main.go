package cmd

import (
	"log"

	"github.com/spf13/cobra"
)

// RootCmd represents the base command when called without any subcommands.
var RootCmd = &cobra.Command{
	Use: "repp",
	Short: `
Repository-based plasmid design. Specify and build plasmids using
their sequence, features, or fragments`,
	Version: "1.0.0",
}

func must(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
