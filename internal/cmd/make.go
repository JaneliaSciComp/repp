package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var (
	backboneHelp = `backbone to insert the fragments into. Can either be an entry 
in one of the dbs or a file on the local filesystem.`

	enzymeHelp = `comma separated list of enzymes to linearize the backbone with.
The backbone must be specified. 'repp ls enzymes' prints a list of
recognized enzymes.`
)

// makeCmd is for finding building a plasmid from its fragments, features, or sequence
var makeCmd = &cobra.Command{
	Use:                        "make",
	Short:                      "Make a plasmid from its expected sequence, features, or fragments",
	SuggestionsMinimumDistance: 3,
	Long: `Find fragments for assembling a plasmid via Gibson Assembly. Build the plasmid
against a list of consituent fragment, feature, or a target sequence.`,
	Aliases: []string{"assemble", "build"},
}

// fragmentsCmd is for piecing together a list of input fragments into a plasmid
var fragmentsCmd = &cobra.Command{
	Use:                        "fragments",
	Short:                      "Build a plasmid from its constituent fragments",
	Run:                        runFragmentsCmd,
	SuggestionsMinimumDistance: 3,
	Long: `Prepare a list of fragments for assembly via Gibson Assembly. Fragments are
checked for existing homology with their neighbors and are prepared for
assembly with PCR.`,
}

// featuresCmd is for building a plasmid from its list of contained features
var featuresCmd = &cobra.Command{
	Use:                        "features \"[feature],...[featureN]\"",
	Short:                      "Find or build a plasmid from its constituent features",
	Run:                        runFeaturesCmd,
	SuggestionsMinimumDistance: 3,
	Example:                    `repp make features "BBa_R0062,BBa_B0034,BBa_C0040,BBa_B0010,BBa_B0012" --backbone pSB1C3 --enzymes "EcoRI,PstI" --dbs igem`,
	Args:                       cobra.MinimumNArgs(1),
}

// sequenceCmd is for assembling a plasmid (single circular sequence) from its target sequence
var sequenceCmd = &cobra.Command{
	Use:                        "sequence",
	Short:                      "Find or build a plasmid from its target sequence",
	Run:                        runSequenceCmd,
	SuggestionsMinimumDistance: 2,
	Long: `Build up a plasmid from its target sequence using a combination of existing and
synthesized fragments.

Solutions have either a minimum fragment count or assembly cost (or both).`,
	Aliases: []string{"seq", "plasmid"},
	Example: `repp make sequence -i "./target_plasmid.fa --dbs addgene`,
}

// set flags
func init() {
	// Flags for specifying the paths to the input file, input fragment files, and output file
	setCommonCommandFlags(fragmentsCmd)
	must(fragmentsCmd.MarkFlagRequired("in"))

	// Flags for specifying the paths to the input file, input fragment files, and output file
	setCommonCommandFlags(featuresCmd)
	must(featuresCmd.MarkFlagRequired("out"))

	// Flags for specifying the paths to the input file, input fragment files, and output file
	setCommonCommandFlags(sequenceCmd)
	must(sequenceCmd.MarkFlagRequired("in"))

	makeCmd.AddCommand(fragmentsCmd)
	makeCmd.AddCommand(featuresCmd)
	makeCmd.AddCommand(sequenceCmd)

	// config is an optional parameter for a settings file (that overrides defaults)
	makeCmd.PersistentFlags().StringP("config", "c", "", "User defined config file that may override all or some default settings")
	makeCmd.PersistentFlags().String("primer3-config", "", "primer3 config folder to be used instead of the default")
	if err := viper.BindPFlag("config", makeCmd.PersistentFlags().Lookup("config")); err != nil {
		log.Fatal(err)
	}

	RootCmd.AddCommand(makeCmd)
}

func setCommonCommandFlags(cmd *cobra.Command) {
	cmd.Flags().StringP("in", "i", "", "input file name (FASTA or Genbank)")
	cmd.Flags().StringP("out", "o", "", "output file name")
	cmd.Flags().StringP("out-fmt", "f", "CSV", "output file format; valid values [JSON, CSV]")
	cmd.Flags().StringP("dbs", "d", "", "list of sequence databases by name")
	cmd.Flags().StringP("enzymes", "e", "", enzymeHelp)
	cmd.Flags().StringP("backbone", "b", "", backboneHelp)
	cmd.Flags().StringP("exclude", "", "", "keywords for excluding fragments")
	cmd.Flags().IntP("identity", "p", 100, "%-identity threshold (see 'blastn -help')")
	cmd.Flags().Bool("ungapped", false, "Ungapped alignment flag")
	cmd.Flags().Int("left-margin", 100, "left margin for matches of the beginning of a circular genome")
	cmd.Flags().StringP("primers-databases", "m", "", "Comma separated list of CSV primers database files")
	cmd.Flags().StringP("primer-prefix", "x", "oS", "Prefix used for labeling new used oligos")
	cmd.Flags().StringP("synth-frags-databases", "s", "", "Comma separated list of CSV synthetic fragments database files")
	cmd.Flags().StringP("synth-frag-prefix", "", "syn", "Prefix used for labeling new found synthetic fragments")
	cmd.Flags().Int("synthetic-frag-factor", 0, "Penalty for synthetic fragments")
	cmd.Flags().IntP("max-kept-solutions", "n", 1, "Top solutions to keep")
}

func runFragmentsCmd(cmd *cobra.Command, args []string) {
	fragmentsInputParams := parseFragmentsAssemblyParams(cmd, args, true)

	if fragmentsInputParams.GetOut() == "" {
		fragmentsInputParams.SetOut(guessOutput(fragmentsInputParams.GetIn(), fragmentsInputParams.GetOutputFormat()))
	}

	syntheticFragmentFactor, err := cmd.Flags().GetInt("synthetic-frag-factor")
	if err != nil {
		log.Printf("Error trying to extract synthetic fragment penalty factor: %v\n", err)
		syntheticFragmentFactor = 0
	}

	config := config.New().SetPrimer3ConfigDir(cmd.Flag("primer3-config").Value.String())
	config.SetSyntheticFragmentFactor(syntheticFragmentFactor)

	repp.AssembleFragments(fragmentsInputParams, config)
}

func runFeaturesCmd(cmd *cobra.Command, args []string) {
	featuresInputParams := parseFeatureAssemblyParams(cmd, args, true)

	if featuresInputParams.GetIn() == "" {
		featuresInputParams.SetIn(combineAllIntoCSV(args))
	}

	syntheticFragmentFactor, err := cmd.Flags().GetInt("synthetic-frag-factor")
	if err != nil {
		log.Printf("Error trying to extract synthetic fragment penalty factor: %v\n", err)
		syntheticFragmentFactor = 0
	}
	maxKeptSolutions, err := cmd.Flags().GetInt("max-kept-solutions")
	if err != nil {
		log.Printf("Error trying to extract synthetic maximum solutions to keep: %v\n", err)
		maxKeptSolutions = 1
	}

	config := config.New().SetPrimer3ConfigDir(cmd.Flag("primer3-config").Value.String())
	config.SetSyntheticFragmentFactor(syntheticFragmentFactor)

	repp.Features(featuresInputParams, maxKeptSolutions, config)
}

func runSequenceCmd(cmd *cobra.Command, args []string) {

	assemblyInputParams := parseSequenceAssemblyParams(cmd, args, true)

	if assemblyInputParams.GetIn() == "" && len(args) > 0 {
		assemblyInputParams.SetIn("input.fa")
		if err := os.WriteFile(assemblyInputParams.GetIn(), []byte(fmt.Sprintf(">target_sequence\n%s", args[0])), 0644); err != nil {
			log.Fatal("Error trying to write target sequence to input.fa", err)
		}
	}

	if assemblyInputParams.GetOut() == "" {
		assemblyInputParams.SetOut(guessOutput(assemblyInputParams.GetIn(), assemblyInputParams.GetOutputFormat()))
	} else {
		assemblyInputParams.SetOut(adjustOutput(assemblyInputParams.GetOut(), assemblyInputParams.GetOutputFormat()))
	}

	syntheticFragmentFactor, err := cmd.Flags().GetInt("synthetic-frag-factor")
	if err != nil {
		log.Printf("Error trying to extract synthetic fragment penalty factor: %v\n", err)
		syntheticFragmentFactor = 0
	}
	maxKeptSolutions, err := cmd.Flags().GetInt("max-kept-solutions")
	if err != nil {
		log.Printf("Error trying to extract synthetic maximum solutions to keep: %v\n", err)
		maxKeptSolutions = 1
	}

	config := config.New().SetPrimer3ConfigDir(cmd.Flag("primer3-config").Value.String())
	config.SetSyntheticFragmentFactor(syntheticFragmentFactor)
	repp.Sequence(assemblyInputParams, maxKeptSolutions, config)
}
