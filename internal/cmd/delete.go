package cmd

import (
	"github.com/spf13/cobra"
)

// deleteCmd is for finding features or enzymes by their name
var deleteCmd = &cobra.Command{
	Use:                        "delete [feature]",
	Short:                      "Delete a feature",
	SuggestionsMinimumDistance: 2,
	Long:                       `Delete a feature, by name, from the embedded feature database.`,
	Aliases:                    []string{"rm", "remove"},
}

// featuresDeleteCmd is for deleting features from the feature db
var featuresDeleteCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Delete a feature from the features database",
	Run:                        featureDB.DeleteCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  repp delete feature \"T7 terminator\"",
	Long: `Delete a feature from the features database by its name.
If no such feature name exists in the database, an error is logged to stderr.`,
}

// set flags
func init() {
	deleteCmd.AddCommand(featuresDeleteCmd)

	RootCmd.AddCommand(deleteCmd)
}
