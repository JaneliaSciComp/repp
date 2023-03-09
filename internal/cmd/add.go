package cmd

import (
	"log"
	"os"
	"path/filepath"
	"strconv"

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
	Short:                      "Import a FASTA sequence database along with its cost",
	Run:                        addDatabase,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nImport a new sequence database so its sequences are available to 'repp make'",
	Example:                    "  repp add database --name addgene --cost 65.0 ./addgene.fa",
	Aliases:                    []string{"db"},
}

// featureAddCmd is for adding a new feature to the features db
var featureAddCmd = &cobra.Command{
	Use:                        "feature [name] [sequence]",
	Short:                      "Add a feature to the features database",
	Run:                        repp.FeaturesAddCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nAdd a feature in the features database so it can be use used in 'repp make features'",
	Example:                    "  repp add feature \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
	Args:                       cobra.ExactArgs(2),
}

// enzymeAddCmd is for adding a new feature to the features db
var enzymeAddCmd = &cobra.Command{
	Use:                        "enzyme [name] [sequence]",
	Short:                      "Add an enzyme to the enzymes database",
	Run:                        repp.EnzymesAddCmd,
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
	databaseAddCmd.Flags().StringP("dir", "d", "", "directory containing FASTA or Genbank file(s) to be used for creating the BLAST database")
	must(databaseAddCmd.MarkFlagRequired("name"))
	must(databaseAddCmd.MarkFlagRequired("cost"))

	addCmd.AddCommand(databaseAddCmd)
	addCmd.AddCommand(featureAddCmd)
	addCmd.AddCommand(enzymeAddCmd)

	RootCmd.AddCommand(addCmd)
}

func addDatabase(cmd *cobra.Command, args []string) {
	dbName := cmd.Flag("name").Value.String()
	costFlag := cmd.Flag("cost").Value.String()
	cost, err := strconv.ParseFloat(costFlag, 64)
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("Cost must be a number", err)
	}
	var seqFiles []string

	seqDir := cmd.Flag("dir").Value.String()
	if seqDir != "" {
		// collect files from sequence dir
		seqDirContent, err := os.ReadDir(seqDir)
		if err != nil {
			log.Fatalf("Error reading directory %s\n", seqDir)
		}
		for _, f := range seqDirContent {
			if !f.IsDir() {
				seqFiles = append(seqFiles, filepath.Join(seqDir, f.Name()))
			}
		}
	}
	// append the rest of files specified as positional args
	seqFiles = append(seqFiles, args...)

	if err = repp.AddDatabase(dbName, seqFiles, cost); err != nil {
		log.Fatal("Error creating database", dbName, err)
	}
}
