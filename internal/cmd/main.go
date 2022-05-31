package cmd

import (
	"github.com/spf13/cobra"
)

// RootCmd represents the base command when called without any subcommands.
var RootCmd = &cobra.Command{
	Use: "repp",
	Short: `REPP
	
Repository-based plasmid design. Specify and build plasmids using
their sequence, features, or fragments`,
	Version: "0.1.0",
}

