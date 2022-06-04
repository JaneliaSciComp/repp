package cmd

import (
	"github.com/Lattice-Automation/repp/internal/repp"

	"github.com/spf13/cobra"
)

// listCmd is for finding features or enzymes by their name.
var listCmd = &cobra.Command{
	Use:                        "list",
	Short:                      "List features or enzymes",
	SuggestionsMinimumDistance: 2,
	Long: `List features or enzymes by name.
If there is no exact match, similar entries are returned`,
	Aliases: []string{"ls"},
}

// databaseListCmd is for reading features (close to the one requested) from the db.
var databaseListCmd = &cobra.Command{
	Use:                        "database [name]",
	Short:                      "List sequence databases",
	Run:                        repp.ListCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  repp list database",
	Long:                       "List all sequence databases and their costs",
	Aliases:                    []string{"db", "dbs", "databases"},
}

// featureListCmd is for reading features (close to the one requested) from the db.
var featureListCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "List features in the features database",
	Run:                        repp.FeaturesReadCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  repp list feature terminator",
	Long: `List features in the features database that are similar to [name].
Writes each feature to the stdout with their name and sequence.

If multiple features contain the feature name sent, each are logged.
Otherwise, all features with names similar to the feature name are writen to stdout`,
}

// enzymeListCmd is for listing out all the available enzymes usable for digesting
// a backbone. Useful for if the user doesn't know which enzymes are available.
var enzymeListCmd = &cobra.Command{
	Use:                        "enzyme [name]",
	Short:                      "List enzymes available for linearizing backbones",
	Run:                        repp.EnzymesReadCmd,
	Example:                    "  repp list enzyme EcoRI",
	SuggestionsMinimumDistance: 2,
	Long: `List out all the enzymes with the same or a similar name as the argument.

'repp list enzyme' without any arguments logs all enzymes available.`,
	Aliases: []string{"enzymes"},
}

// fragmentListCmd is for finding a fragment by its name
var fragmentListCmd = &cobra.Command{
	Use:                        "fragment [name]",
	Short:                      "List a fragment in the databases",
	Example:                    "  repp list fragment pSB1C3 --dbs igem",
	Run:                        repp.FragmentListCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       `List fragments with a passed name in the specified databases`,
}

// sequenceListCmd is for finding a sequence in the dbs
var sequenceListCmd = &cobra.Command{
	Use:                        "sequence [seq]",
	Short:                      "List a sequence in the databases",
	Run:                        repp.SequenceListCmd,
	Example:                    "  repp list sequence GTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGAC --dbs igem",
	SuggestionsMinimumDistance: 2,
	Long:                       `List a sequence's BLAST matches among databases.`,
	Aliases:                    []string{"seq"},
}

// set flags
func init() {
	fragmentListCmd.Flags().StringP("dbs", "d", "", "comma separated list of sequence databases")

	sequenceListCmd.Flags().StringP("dbs", "d", "", "comma separated list of sequence databases")
	sequenceListCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	sequenceListCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	listCmd.AddCommand(databaseListCmd)
	listCmd.AddCommand(featureListCmd)
	listCmd.AddCommand(enzymeListCmd)
	listCmd.AddCommand(fragmentListCmd)
	listCmd.AddCommand(sequenceListCmd)

	RootCmd.AddCommand(listCmd)
}
