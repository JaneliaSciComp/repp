package repp

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/spf13/cobra"
)

// Annotate is for annotating a plasmid sequence given the features in the feature database.
// If an output path is provided, the annotated plasmid is writen to that file. Otherwise,
// the feature matches are written to stdout.
func Annotate(cmd *cobra.Command, args []string) {
	output, _ := cmd.Flags().GetString("out")

	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = 96 // might be something other than `repp plasmid`
	}

	p := inputParser{}
	filters, err := cmd.Flags().GetString("exclude")
	excludeFilters := []string{}
	if err == nil {
		excludeFilters = p.getFilters(filters)
	}

	name := ""
	query := "" // the plasmid sequence that we're querying
	if len(args) > 0 {
		query = args[0]
	} else {
		in, err := cmd.Flags().GetString("in")
		if in == "" || err != nil {
			if helperr := cmd.Help(); helperr != nil {
				log.Fatal(helperr)
			}
			stderr.Fatalln("must pass a file with a plasmid sequence or the plasmid sequence as an argument.")
		}

		frags, err := read(in, false)
		if err != nil {
			stderr.Fatalln(err)
		}
		name = frags[0].ID
		query = frags[0].Seq
	}

	toCull, _ := cmd.Flags().GetBool("cull")
	namesOnly, _ := cmd.Flags().GetBool("names")

	dbflag, err := cmd.Flags().GetString("dbs")
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		stderr.Fatalf("failed to parse building fragments: %v", err)
	}

	m, err := newManifest()
	if err != nil {
		stderr.Fatalf("failed to get DB manifest: %v", err)
	}
	dbs, err := p.parseDBs(m, dbflag)
	if err != nil {
		stderr.Fatalf("failed to find any fragment databases: %v", err)
	}

	annotate(name, query, output, identity, dbs, excludeFilters, toCull, namesOnly)
}

// annotate is for executing blast against the query sequence.
func annotate(name, seq, output string, identity int, dbs []DB, filters []string, toCull, namesOnly bool) {
	handleErr := func(err error) {
		if err != nil {
			stderr.Fatalln(err)
		}
	}

	in, err := ioutil.TempFile("", "annotate-in-*")
	handleErr(err)
	defer os.Remove(in.Name())

	out, err := ioutil.TempFile("", "annotate-out-*")
	handleErr(err)
	defer os.Remove(out.Name())

	// create a subject file with all the blast features
	featureKV := NewFeatureDB()
	featIndex := 0
	var featureSubjects strings.Builder
	indexToFeature := make(map[int]string)
	for feat, featSeq := range featureKV.contents {
		indexToFeature[featIndex] = feat
		featureSubjects.WriteString(fmt.Sprintf(">%d\n%s\n", featIndex, featSeq))
		featIndex++
	}
	subjectFile, err := ioutil.TempFile("", "features-*")
	handleErr(err)
	defer os.Remove(subjectFile.Name())

	_, err = subjectFile.WriteString(featureSubjects.String())
	handleErr(err)

	b := &blastExec{
		in:       in,
		out:      out,
		name:     name,
		subject:  subjectFile.Name(),
		seq:      seq,
		identity: identity,
		circular: true,
	}

	var features []match
	if len(dbs) < 1 {
		// if the user selected another db, don't use the internal one
		handleErr(b.input())
		handleErr(b.runAgainst())
		features, err = b.parse(filters)
		handleErr(err)

		// get rid of features that start past the zero index, wrap that those that go around it
		// get rid of features matches that aren't 100% of the feature in the feature database
		var cleanedFeatures []match
		for _, f := range features {
			if f.queryStart >= len(seq) {
				continue
			}

			featureIndex, _ := strconv.Atoi(f.entry)
			f.entry = indexToFeature[featureIndex]
			if len(f.seq) < len(featureKV.contents[f.entry]) {
				continue
			}

			f.queryEnd %= len(seq)
			if f.queryEnd == 0 {
				f.queryEnd = len(seq)
			}

			cleanedFeatures = append(cleanedFeatures, f)
		}
		features = cleanedFeatures
	} else {
		features, err = blast(name, seq, false, dbs, filters, identity, blastWriter())
		handleErr(err)
	}

	if len(features) < 1 {
		stderr.Fatal("no features found")
	}

	sortMatches(features)
	if toCull {
		features = cull(features, len(seq), 5, 1)
		// for _, f := range features {
		// 	f.log()
		// }
	}

	if namesOnly {
		featuresNames := []string{}
		for _, feature := range features {
			dir := ""
			if !feature.forward {
				dir += ":rev"
			}
			featuresNames = append(featuresNames, feature.entry+dir)
		}
		fmt.Println(strings.Join(featuresNames, ", "))
	} else if output != "" {
		writeGenbank(output, name, seq, []*Frag{}, features)
	} else {
		tw := tabwriter.NewWriter(os.Stdout, 0, 4, 3, ' ', 0)
		fmt.Fprintf(tw, "\nfeatures (%d)\tstart\tend\tdirection\t\n", len(features))
		for _, feat := range features {
			dir := "FWD"
			if !feat.forward {
				dir = "REV"
			}
			fmt.Fprintf(tw, "%s\t%d\t%d\t%s\t\n", feat.entry, feat.queryStart+1, feat.queryEnd+1, dir)
		}
		tw.Flush()
	}
}
