package cmd

import (
	"log"
	"strings"

	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

// deleteCmd is for finding features or enzymes by their name
var deleteCmd = &cobra.Command{
	Use:                        "delete [database,feature]",
	Short:                      "Delete a feature",
	SuggestionsMinimumDistance: 2,
	Long:                       `Delete a feature, by name, from the embedded feature database.`,
	Aliases:                    []string{"rm", "remove"},
}

// databaseDeleteCmd is for deleting features from the feature db
var databaseDeleteCmd = &cobra.Command{
	Use:                        "database [name]",
	Short:                      "Delete a sequence database",
	Run:                        runDatabaseDeleteCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  repp delete database \"igem\"",
	Aliases:                    []string{"db"},
	Args:                       cobra.ExactArgs(1),
}

// featuresDeleteCmd is for deleting features from the feature db
var featuresDeleteCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Delete a feature from the features database",
	Run:                        runFeaturesDeleteCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  repp delete feature \"T7 terminator\"",
	Long: `Delete a feature from the features database by its name.
If no such feature name exists in the database, an error is logged to stderr.`,
	Args: cobra.ExactArgs(1),
}

// set flags
func init() {
	deleteCmd.AddCommand(databaseDeleteCmd)
	deleteCmd.AddCommand(featuresDeleteCmd)

	RootCmd.AddCommand(deleteCmd)
}

func runDatabaseDeleteCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("No database was specified")
	}
	db := args[0]

	repp.DeleteDatabase(db)
}

func runFeaturesDeleteCmd(cmd *cobra.Command, args []string) {
	var name string

	if len(args) < 1 {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("\nno features name passed.")
	} else if len(args) == 1 {
		name = args[0]
	} else {
		name = strings.Join(args, " ")
	}

	repp.DeleteFeature(name)
}
