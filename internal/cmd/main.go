package cmd

import (
	"log"

	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

// RootCmd represents the base command when called without any subcommands.
var RootCmd = &cobra.Command{
	Use:   "repp",
	Short: `repository-based plasmid design. Build cost-efficient plasmids`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		if cmd.Flag("verbose").Value.String() == "true" {
			repp.SetVerboseLogging()
		}
	},
	Version: "1.0.0",
}

func init() {
	RootCmd.PersistentFlags().BoolP("verbose", "v", false, "write DEBUG logs")
}

func must(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
