package cmd

import (
	"log"
	"strings"

	"github.com/Lattice-Automation/repp/internal/repp"
	"github.com/spf13/cobra"
	"golang.org/x/exp/slices"
)

// ParseFeatureAssemblyParams - parse feature specific flags from the command line
func ParseFeatureAssemblyParams(cmd *cobra.Command, args []string, strict bool) *repp.AssemblyParams {
	params := &repp.AssemblyParams{}

	extractCommonParams(cmd, args, params)
	// extract filters
	params.Filters = extractExcludedValues(cmd)

	return params
}

// ParseSequenceAssemblyParams - parse sequence specific flags from the command line
func ParseSequenceAssemblyParams(cmd *cobra.Command, args []string, strict bool) *repp.AssemblyParams {
	params := &repp.AssemblyParams{}

	extractCommonParams(cmd, args, params)
	// extract filters
	params.Filters = extractExcludedValues(cmd)

	return params
}

// ParseFragmentsAssemblyParams - parse sequence specific flags from the command line
func ParseFragmentsAssemblyParams(cmd *cobra.Command, args []string, strict bool) *repp.AssemblyParams {
	params := &repp.AssemblyParams{}

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

func extractCommonParams(cmd *cobra.Command, args []string, params *repp.AssemblyParams) {
	inputFname, err := cmd.Flags().GetString("in")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse sequence input arg: %v", err)
	}
	params.In = inputFname

	outputFName, err := cmd.Flags().GetString("out")

	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse output arg: %v", err)
	}
	params.Out = outputFName

	// get identity for blastn searching
	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = 100 // might be something other than `repp plasmid`
	}
	params.Identity = identity

	dbNames, err := cmd.Flags().GetString("dbs")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatalf("failed to parse dbs arg: %v", err)
	}
	params.DbNames = splitStringOn(dbNames, []rune{' ', ','})

	// check if user asked for a specific backbone, confirm it exists in one of the dbs
	backbone, _ := cmd.Flags().GetString("backbone")
	params.BackboneName = backbone

	// check if they also specified an enzyme
	enzymes, _ := cmd.Flags().GetString("enzymeList")
	params.EnzymeNames = splitStringOn(enzymes, []rune{' ', ','})
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
