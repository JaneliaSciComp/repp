package cmd

import (
	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

var (
	featureDB = repp.NewFeatureDB()

	enzymeDB = repp.NewEnzymeDB()
)

// RootCmd represents the base command when called without any subcommands.
var RootCmd = &cobra.Command{
	Use: "repp",
	Short: `REPP
	
Repository-based plasmid design. Specify and build plasmids using
their sequence, features, or fragments`,
	Version: "0.1.0",
}

