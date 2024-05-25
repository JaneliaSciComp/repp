package repp

import (
	"fmt"
	"os"
	"strconv"
	"strings"
	"text/tabwriter"
)

// Annotate is for annotating a plasmid sequence given the features in the feature database.
// If an output path is provided, the annotated plasmid is writen to that file. Otherwise,
// the feature matches are written to stdout.
func Annotate(inputName, inputQuery string,
	identity int,
	namesOnly, toCull bool,
	dbNames, filters []string,
	output string) {
	var name, query string

	if inputQuery == "" {
		if inputName == "" {
			rlog.Fatal("must pass a file with a plasmid sequence or the plasmid sequence as an argument.")
		} else {
			frags, err := read(inputName, false, false)
			if err != nil {
				rlog.Fatal(err)
			}
			name = frags[0].ID
			query = frags[0].Seq
		}
	} else {
		name = inputName
		query = inputQuery
	}

	dbs, err := getRegisteredDBs(dbNames)
	if err != nil {
		rlog.Fatal("failed to find any fragment databases: %v", err)
	}

	annotate(name, query, output, identity, dbs, filters, toCull, namesOnly)
}

// annotate is for executing blast against the query sequence.
func annotate(name, seq, output string, identity int, dbs []DB, filters []string, toCull, namesOnly bool) {
	handleErr := func(err error) {
		if err != nil {
			rlog.Fatal(err)
		}
	}

	in, err := os.CreateTemp("", "annotate-in-*")
	handleErr(err)

	out, err := os.CreateTemp("", "annotate-out-*")
	handleErr(err)

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
	subjectFile, err := os.CreateTemp("", "features-*")
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
		ungapped: false,
	}
	defer b.close()

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
		features, err = blast(name, seq, false, 0, dbs, filters, identity, false)
		handleErr(err)
	}

	if len(features) < 1 {
		rlog.Fatal("no features found")
	}

	sortMatches(features)
	if toCull {
		features = cull(features, 5, 1)
	}

	if namesOnly {
		featuresNames := []string{}
		for _, feature := range features {
			dir := ""
			if feature.isRevCompMatch() {
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
			if feat.isRevCompMatch() {
				dir = "REV"
			}
			fmt.Fprintf(tw, "%s\t%d\t%d\t%s\t\n", feat.entry, feat.queryStart+1, feat.queryEnd+1, dir)
		}
		tw.Flush()
	}
}
