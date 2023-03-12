package cmd

import (
	"log"
	"path/filepath"
	"strings"

	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
	"golang.org/x/exp/slices"
)

// parseFeatureAssemblyParams - parse feature specific flags from the command line
func parseFeatureAssemblyParams(cmd *cobra.Command, args []string, strict bool) repp.AssemblyParams {
	params := repp.MkAssemblyParams()

	extractCommonParams(cmd, args, params)
	// extract filters
	params.SetFilters(extractExcludedValues(cmd))

	return params
}

// parseSequenceAssemblyParams - parse sequence specific flags from the command line
func parseSequenceAssemblyParams(cmd *cobra.Command, args []string, strict bool) repp.AssemblyParams {
	params := repp.MkAssemblyParams()

	extractCommonParams(cmd, args, params)
	// extract filters
	params.SetFilters(extractExcludedValues(cmd))

	return params
}

// parseFragmentsAssemblyParams - parse sequence specific flags from the command line
func parseFragmentsAssemblyParams(cmd *cobra.Command, args []string, strict bool) repp.AssemblyParams {
	params := repp.MkAssemblyParams()

	extractCommonParams(cmd, args, params)

	return params
}

func extractExcludedValues(cmd *cobra.Command) []string {
	excluded, err := cmd.Flags().GetString("exclude")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse excluded arg: %v", err)
	}

	return splitStringOn(strings.ToUpper(excluded), []rune{' ', ','})
}

func extractIdentity(cmd *cobra.Command, defaultValue int) int {
	// get identity for blastn searching
	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = defaultValue // might be something other than `repp plasmid`
	}
	return identity
}

func extractDbNames(cmd *cobra.Command) []string {
	dbNames, err := cmd.Flags().GetString("dbs")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse dbs arg: %v", err)
	}
	return splitStringOn(dbNames, []rune{' ', ','})
}

func extractCommonParams(cmd *cobra.Command, args []string, params repp.AssemblyParams) {
	inputFname, err := cmd.Flags().GetString("in")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse sequence input arg: %v", err)
	}
	params.SetIn(inputFname)

	outputFName, err := cmd.Flags().GetString("out")

	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse output arg: %v", err)
	}
	params.SetOut(outputFName)

	// get identity for blastn searching
	params.SetIdentity(extractIdentity(cmd, 100))

	params.SetDbNames(extractDbNames(cmd))

	// check if user asked for a specific backbone, confirm it exists in one of the dbs
	backboneName, _ := cmd.Flags().GetString("backbone")
	params.SetBackboneName(backboneName)

	// check if they also specified an enzyme
	enzymeNames, _ := cmd.Flags().GetString("enzymeList")
	params.SetEnzymeNames(splitStringOn(enzymeNames, []rune{' ', ','}))
}

// guessOutput gets an outpath path from an input path (if no output path is
// specified). It uses the same name as the input path to create an output.
func guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".output.json"
}

func splitStringOn(s string, separators []rune) []string {
	splitFunc := func(c rune) bool {
		return slices.Contains(separators, c)
	}

	return strings.FieldsFunc(s, splitFunc)
}

// combineAllIntoCSV turns the arguments into a
// CSV list of values in a single string.
func combineAllIntoCSV(args []string) string {
	commaSeparated := false
	for _, a := range args {
		if strings.Contains(a, ",") {
			commaSeparated = true
		}
	}

	// if it's a features command, concatenate the arguments in case they're feature names
	// with 'repp features' the arguments are the feature names to use
	if commaSeparated {
		spacedSeparated := strings.Join(args, " ")
		splitByComma := strings.Split(spacedSeparated, ",")
		trimmedSeqs := []string{}
		for _, entry := range splitByComma {
			trimmedSeqs = append(trimmedSeqs, strings.TrimSpace(entry))
		}
		return strings.Join(trimmedSeqs, ",")
	}

	return strings.Join(args, " ")
}
