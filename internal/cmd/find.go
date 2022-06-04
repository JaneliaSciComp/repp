package cmd

import (
	"Lattice-Automation/repp/internal/repp"

	"github.com/spf13/cobra"
)

// findCmd is for finding features or enzymes by their name.
var findCmd = &cobra.Command{
	Use:                        "find",
	Short:                      "Find features or enzymes",
	SuggestionsMinimumDistance: 2,
	Long: `Find features or enzymes by name.
If there is no exact match, similar entries are returned`,
	Aliases: []string{"ls", "list", "get"},
}

// featureFindCmd is for reading features (close to the one requested) from the db.
var featureFindCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Find features in the features database",
	Run:                        featureDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  repp find feature terminator",
	Long: `Find features in the features database that are similar to [name].
Writes each feature to the stdout with their name and sequence.

If multiple features contain the feature name sent, each are logged.
Otherwise, all features with names similar to the feature name are writen to stdout`,
}

// enzymeFindCmd is for listing out all the available enzymes usable for digesting
// a backbone. Useful for if the user doesn't know which enzymes are available.
var enzymeFindCmd = &cobra.Command{
	Use:                        "enzyme [name]",
	Short:                      "Find enzymes available for linearizing backbones",
	Run:                        enzymeDB.ReadCmd,
	Example:                    "  repp find enzyme EcoRI",
	SuggestionsMinimumDistance: 2,
	Long: `List out all the enzymes with the same or a similar a similar name as the argument.

'repp find enzyme' without any arguments logs all enzymes available.`,
	Aliases: []string{"enzymes"},
}

// fragmentFindCmd is for finding a fragment by its name
var fragmentFindCmd = &cobra.Command{
	Use:                        "fragment [name]",
	Short:                      "Find a fragment in the databases",
	Example:                    "  repp find fragment pSB1C3 --igem",
	Run:                        repp.FragmentFindCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       `Find a fragment with a given name in the databases requested.`,
}

// sequenceFindCmd is for finding a sequence in the dbs
var sequenceFindCmd = &cobra.Command{
	Use:                        "sequence [seq]",
	Short:                      "Find a sequence in the databases",
	Run:                        repp.SequenceFindCmd,
	Example:                    "  repp find sequence GTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGAC --igem",
	SuggestionsMinimumDistance: 2,
	Long:                       `Find a sequence's BLAST matches among databases.`,
	Aliases:                    []string{"seq"},
}

// set flags
func init() {
	fragmentFindCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases")
	fragmentFindCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	fragmentFindCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	fragmentFindCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU respository")

	sequenceFindCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases")
	sequenceFindCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	sequenceFindCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	sequenceFindCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU respository")
	sequenceFindCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	sequenceFindCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	findCmd.AddCommand(featureFindCmd)
	findCmd.AddCommand(enzymeFindCmd)
	findCmd.AddCommand(fragmentFindCmd)
	findCmd.AddCommand(sequenceFindCmd)

	RootCmd.AddCommand(findCmd)
}
