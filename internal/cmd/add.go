package cmd

import (
	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

// addCmd is for piecing together a list of input fragments into a plasmid
// and preparing the fragments to make into that plasmid
var addCmd = &cobra.Command{
	Use:                        "add [database,feature,enzyme]",
	Short:                      "Add a sequence database, feature, or enzyme",
	SuggestionsMinimumDistance: 1,
	Long: `Create/update a feature or enzyme with its name and sequence/recognition-site.
Set features can be passed to the 'repp build features' command and enzymes can
be passed to the --enzyme flag`,
	Aliases: []string{"set", "update"},
}

// databaseAddCmd is for adding a new sequence db
var databaseAddCmd = &cobra.Command{
	Use:                        "database [path] [cost-per-seq]",
	Short:                      "Import a FASTA sequence database along with its cost",
	Run:                        repp.AddCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nImport a new sequence database so its sequences are available to 'repp make'",
	Example:                    "  repp add database ./addgene.fa 65.0",
	Aliases: []string{"db"},
}

// featureAddCmd is for adding a new feature to the features db
var featureAddCmd = &cobra.Command{
	Use:                        "feature [name] [sequence]",
	Short:                      "Add a feature to the features database",
	Run:                        repp.FeaturesAddCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nAdd a feature in the features database so it can be use used in 'repp make features'",
	Example:                    "  repp add feature \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
}

// enzymeAddCmd is for adding a new feature to the features db
var enzymeAddCmd = &cobra.Command{
	Use:                        "enzyme [name] [sequence]",
	Short:                      "Add an enzyme to the enzymes database",
	Run:                        repp.EnzymesAddCmd,
	SuggestionsMinimumDistance: 2,
	Long: `Add an enzyme in the enzymes database so it can be used to linearize backbones.
Enzymes are passed to the build command, by name, with the --enzyme flag.

Valid recognition sequences have both a cut site in the template sequence: "^" and
a cut site in the complement sequence: "_". Use 'repp ls enzyme' for examples`,
	Example: "  repp add enzyme BbvCI CC^TCA_GC",
}

func init() {
	addCmd.AddCommand(databaseAddCmd)
	addCmd.AddCommand(featureAddCmd)
	addCmd.AddCommand(enzymeAddCmd)

	RootCmd.AddCommand(addCmd)
}
