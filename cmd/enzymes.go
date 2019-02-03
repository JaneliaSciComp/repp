package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// enzymesCmd is for listing out all the available enzymes usable for digesting
// a backbone. Useful for if the user doesn't know which enzymes are available
var enzymesCmd = &cobra.Command{
	Use:                        "enzymes",
	Short:                      "List enzymes available to linearize a backbone",
	Run:                        defrag.Enzymes,
	SuggestionsMinimumDistance: 4,
	Long: `Lists out all the enzymes in defrag by name along with their recognition sequence.

<Name>: <Recognition sequence>`,
}

func init() {
	rootCmd.AddCommand(enzymesCmd)
}
