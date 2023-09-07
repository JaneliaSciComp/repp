package cmd

import (
	"log"
	"strings"

	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

// addCmd is for piecing together a list of input fragments into a plasmid
// and preparing the fragments to make into that plasmid
var addCmd = &cobra.Command{
	Use:                        "add",
	Short:                      "Add a sequence database, feature, or enzyme",
	SuggestionsMinimumDistance: 1,
	Long: `Create/update a feature or enzyme with its name and sequence/recognition-site.
Set features can be passed to the 'repp build features' command and enzymes can
be passed to the --enzyme flag`,
	Aliases: []string{"set"},
}

// databaseAddCmd is for adding a new sequence db
var databaseAddCmd = &cobra.Command{
	Use:                        "database",
	Short:                      "Import a FASTA sequence database along with its cost.",
	Run:                        runDatabaseAddCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nImport a new sequence database so its sequences are available to 'repp make'.",
	Example:                    "  repp add database --name addgene --cost 65.0 ./addgene.fa",
	Aliases:                    []string{"db"},
}

// featureAddCmd is for adding a new feature to the features db
var featureAddCmd = &cobra.Command{
	Use:                        "feature [name] [sequence]",
	Short:                      "Add a feature to the features database",
	Run:                        runFeaturesAddCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nAdd a feature in the features database so it can be use used in 'repp make features'",
	Example:                    "  repp add feature \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
	Args:                       cobra.ExactArgs(2),
}

// enzymeAddCmd is for adding a new feature to the features db
var enzymeAddCmd = &cobra.Command{
	Use:                        "enzyme [name] [sequence]",
	Short:                      "Add an enzyme to the enzymes database",
	Run:                        runEnzymesAddCmd,
	SuggestionsMinimumDistance: 2,
	Long: `Add an enzyme in the enzymes database so it can be used to linearize backbones.
See: 'repp make sequence --help' for usage of enzymes.

Valid recognition sequences have both a cut site in the template sequence: "^" and
a cut site in the complement sequence: "_". Use 'repp ls enzyme' for examples`,
	Example: "  repp add enzyme BbvCI CC^TCA_GC",
	Args:    cobra.ExactArgs(2),
}

func init() {
	databaseAddCmd.Flags().StringP("name", "n", "", "database name")
	databaseAddCmd.Flags().Float64P("cost", "c", 0.0, "the cost per plasmid procurement (eg order + shipping fee)")
	databaseAddCmd.Flags().Bool("prefixSeqIDs", true, "Prefix sequence IDs with filename")

	must(databaseAddCmd.MarkFlagRequired("name"))
	must(databaseAddCmd.MarkFlagRequired("cost"))

	addCmd.AddCommand(databaseAddCmd)
	addCmd.AddCommand(featureAddCmd)
	addCmd.AddCommand(enzymeAddCmd)

	RootCmd.AddCommand(addCmd)
}

func runDatabaseAddCmd(cmd *cobra.Command, args []string) {
	dbName, err := cmd.Flags().GetString("name")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("Database name must be a string", err)
	}
	cost, err := cmd.Flags().GetFloat64("cost")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("Cost must be a number", err)
	}
	prefixSeqIDs, err := cmd.Flags().GetBool("prefixSeqIDs")
	if err != nil {
		log.Print("Error encountered reading prefiSeqIDs flag", err)
		prefixSeqIDs = false
	}

	seqFiles, err := repp.CollectFiles(args)
	if err != nil {
		log.Fatalf("Errors encountered collection sequence files from %v: %v", args, err)
	}

	if err = repp.AddDatabase(dbName, seqFiles, cost, prefixSeqIDs); err != nil {
		log.Fatal("Error creating database", dbName, err)
	}
}

func runFeaturesAddCmd(cmd *cobra.Command, args []string) {
	var name, seq string

	if len(args) < 2 {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("Add features must have exactly 2 arguments")
		return
	} else if len(args) == 2 {
		name = args[0]
		seq = args[1]
	} else {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]

	}

	repp.AddFeatures(name, seq)
}

func runEnzymesAddCmd(cmd *cobra.Command, args []string) {
	var name, seq string

	if len(args) < 2 {

	} else if len(args) == 2 {
		name = args[0]
		seq = args[1]
	} else {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]
	}

	repp.AddEnzymes(name, seq)
}
