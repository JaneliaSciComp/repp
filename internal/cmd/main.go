package cmd

import (
	_ "embed"
	"fmt"
	"log"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

const releaseNumber = "1.1.0"

//go:embed commit.txt
var commit string

// RootCmd represents the base command when called without any subcommands.
var RootCmd = &cobra.Command{
	Use:   "repp",
	Short: `repository-based plasmid design. Build cost-efficient plasmids`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		if cmd.Flag("verbose").Value.String() == "true" {
			repp.SetVerboseLogging()
		}
		reppDataDir := cmd.Flag("repp-data-dir").Value.String()

		config.Setup(reppDataDir)
	},
	Version: fmt.Sprintf("%s (%.11s)", releaseNumber, commit),
}

func init() {
	RootCmd.PersistentFlags().BoolP("verbose", "v", false, "write DEBUG logs")
	RootCmd.PersistentFlags().String("repp-data-dir", "", "Default REPP data directory")
}

func must(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
