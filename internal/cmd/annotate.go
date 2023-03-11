package cmd

import (
	"log"

	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
)

// annotateCmd is for finding features or enzymes by their name.
var annotateCmd = &cobra.Command{
	Use:                        "annotate [seq]",
	Run:                        runAnnotateCmd,
	Short:                      "Annotate a plasmid using features",
	SuggestionsMinimumDistance: 3,
	Long: `Accepts a sequence file as input and runs alignment against the
embedded feature database. Each alignment feature is included as
a feature in the output: a Genbank file. Individual databases
can be selected, in which case the entries in the database will
be used in the alignment _rather_ than the feature database.

The feature database and the default 96% identity are based on
information from [SnapGene](https://www.snapgene.com/resources/plasmid-files/)`,
}

// set flags
func init() {
	annotateCmd.Flags().StringP("in", "i", "", "input file name")
	annotateCmd.Flags().StringP("out", "o", "", "output file name")
	annotateCmd.Flags().StringP("exclude", "x", "", "keywords for excluding features")
	annotateCmd.Flags().StringP("dbs", "d", "", "comma separated list sequence databases to consider as features")
	annotateCmd.Flags().IntP("identity", "p", 96, "match %-identity threshold (see 'blastn -help')")
	annotateCmd.Flags().BoolP("cull", "c", true, "remove features enclosed in others")
	annotateCmd.Flags().BoolP("names", "n", false, "log feature names to the console")

	RootCmd.AddCommand(annotateCmd)
}

func runAnnotateCmd(cmd *cobra.Command, args []string) {
	var name string
	var query string

	name, _ = cmd.Flags().GetString("in")

	output, _ := cmd.Flags().GetString("out")

	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = 96 // might be something other than `repp plasmid`
	}

	filters := extractExcludedValues(cmd)

	if len(args) > 0 {
		query = args[0]
	}

	if name == "" && query == "" {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("must pass a file with a plasmid sequence or the plasmid sequence as an argument.")
	}

	namesOnly, _ := cmd.Flags().GetBool("names")
	toCull, _ := cmd.Flags().GetBool("cull")

	dbNamesValue, err := cmd.Flags().GetString("dbs")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse dbs arg: %v", err)
	}
	dbNames := splitStringOn(dbNamesValue, []rune{' ', ','})

	repp.Annotate(
		name,
		query,
		identity,
		namesOnly,
		toCull,
		dbNames,
		filters,
		output)
}
