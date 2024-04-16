package repp

import (
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
	"time"

	"github.com/Lattice-Automation/repp/internal/config"
)

type featureMatch struct {
	featureIndex int
	match        match
}

// Features assembles a plasmid with all the Features requested with the 'repp Features [feature ...]' command
// repp assemble Features p10 promoter, mEGFP, T7 terminator
func Features(assemblyParams AssemblyParams, maxSolutions int, conf *config.Config) [][]*Frag {
	start := time.Now()

	// get registered blast databases
	dbs, err := assemblyParams.getDBs()
	if err != nil {
		// error getting the DBs
		rlog.Fatal(err)
	}
	// get registered enzymes
	enzymes, err := assemblyParams.getEnzymes()
	if err != nil {
		// error getting the enzymes
		rlog.Fatal(err)
	}
	// prepare backbone if needed
	backboneFrag, backboneMeta, err := prepareBackbone(assemblyParams.GetBackboneName(), enzymes, dbs)
	if err != nil {
		// error getting the backbone
		rlog.Fatal(err)
	}

	// turn feature names into sequences
	insertFeats, bbFeat := queryFeatures(
		assemblyParams.GetIn(),
		backboneFrag,
		dbs,
	)
	feats := insertFeats
	if len(bbFeat) > 0 {
		feats = append(feats, bbFeat)
	}

	// find matches in the databases
	featureMatches := blastFeatures(
		assemblyParams.GetFilters(),
		assemblyParams.GetIdentity(),
		dbs,
		feats,
		conf,
	)
	if len(featureMatches) == 0 {
		featNames := []string{}
		for _, feat := range insertFeats {
			featNames = append(featNames, feat[0])
		}
		rlog.Fatal("failed to find fragments with specified features", "features", featNames)
	}

	// build assemblies containing the matched fragments
	target, solutions := featureSolutions(
		feats,
		featureMatches,
		assemblyParams.GetIdentity(),
		dbs,
		maxSolutions,
		conf)

	// write the output file
	insertLength := 0
	for _, f := range insertFeats {
		insertLength += len(f[1])
	}

	// do not use the oligos manifest
	primersDB := readOligos(assemblyParams.GetPrimersDBLocations(), primerIDPrefix, false)
	synthFragsDB := readOligos(assemblyParams.GetSynthFragsDBLocations(), synthFragIDPrefix, true)

	if _, err := writeResult(
		assemblyParams.GetOut(),
		assemblyParams.GetOutputFormat(),
		assemblyParams.GetIn(),
		target,
		solutions,
		primersDB,
		synthFragsDB,
		time.Since(start).Seconds(),
		backboneMeta,
		conf,
	); err != nil {
		rlog.Fatal(err)
	}

	return solutions
}

// queryFeatures takes the list of feature names and finds them in the available databases
func queryFeatures(
	featuresInput string,
	backbone *Frag,
	dbs []DB) ([][]string, []string) {
	var insertFeats [][]string // slice of tuples [feature name, feature sequence]
	if readFeatures, err := read(featuresInput, true, false); err == nil {
		// see if the features are in a file (multi-FASTA or features in a Genbank)
		seenFeatures := make(map[string]string) // map feature name to sequence
		for _, f := range readFeatures {
			if seq := seenFeatures[f.ID]; seq != f.Seq {
				rlog.Fatal("failed to parse features, %s has two different sequences:\n\t%s\n\t%s\n", f.ID, f.Seq, seq)
			}
			insertFeats = append(insertFeats, []string{f.ID, f.Seq})
		}
	} else {
		// if the features weren't in a file, try and find each in the features database
		// or one of the databases passed as a source of building fragments
		var featureNames []string
		if strings.Contains(featuresInput, ",") {
			featureNames = strings.Split(featuresInput, ",") // comma separated
		} else {
			featureNames = strings.Fields(featuresInput) // spaces
		}

		if len(featureNames) < 1 {
			rlog.Fatal("no features chosen. see 'repp make features --help'")
		}

		featureDB := NewFeatureDB()
		for _, f := range featureNames {
			fwd := true
			if strings.Contains(f, ":") {
				ns := strings.Split(f, ":")
				f = ns[0]
				fwd = !strings.Contains(strings.ToLower(ns[1]), "rev")
			}

			if seq, contained := featureDB.contents[f]; contained {
				if !fwd {
					f = f + ":REV"
					seq = reverseComplement(seq)
				}
				insertFeats = append(insertFeats, []string{f, seq})
			} else if dbFrag, err := queryDatabases(f, dbs); err == nil {
				f = strings.Replace(f, ":", "|", -1)
				if !fwd {
					dbFrag.Seq = reverseComplement(dbFrag.Seq)
				}
				insertFeats = append(insertFeats, []string{f, dbFrag.Seq})
			} else {
				rlog.Fatalf(
					"failed to find '%s' among the features in (%s) or any db: %s",
					f,
					config.FeatureDB,
					strings.Join(dbNames(dbs), ","),
				)
			}
		}
	}

	// add in the backbone as a "feature"
	bbFeat := []string{}
	if backbone != nil && backbone.ID != "" {
		bbFeat = []string{backbone.ID, backbone.Seq}
	}

	return insertFeats, bbFeat
}

// blastFeatures returns matches between the target features and entries in the databases with those features
func blastFeatures(
	filters []string,
	identity int,
	dbs []DB,
	feats [][]string,
	conf *config.Config) map[string][]featureMatch {
	featureMatches := make(map[string][]featureMatch) // a map from from each entry (by id) to its list of matched features
	for i, target := range feats {
		targetFeature := target[1]
		matches, err := blast(
			target[0],
			targetFeature,
			false,
			0,
			dbs,
			filters,
			identity)
		if err != nil {
			rlog.Fatal(err)
		}

		for _, m := range matches {
			// needs to be at least identity % as long as the queried feature
			mLen := float64(m.subjectEnd - m.subjectStart + 1)
			pIdent := mLen / float64(len(targetFeature))
			pIdentTarget := float64(identity) / 100.0
			if pIdent < pIdentTarget {
				continue
			}

			m.queryStart = i
			m.queryEnd = i
			m.uniqueID = m.entry + strconv.Itoa(m.subjectStart)

			if _, exists := featureMatches[m.entry]; !exists {
				featureMatches[m.entry] = []featureMatch{{featureIndex: i, match: m}}
			} else {
				featureMatches[m.entry] = append(featureMatches[m.entry], featureMatch{featureIndex: i, match: m})
			}
		}
	}

	return featureMatches
}

// featureSolutions creates and fills the assemblies using the matched fragments
func featureSolutions(
	feats [][]string,
	featureMatches map[string][]featureMatch,
	identity int,
	dbs []DB,
	keepNSolutions int,
	conf *config.Config) (string, [][]*Frag) {
	// merge matches into one another if they can combine to cover a range
	extendedMatches := extendMatches(feats, featureMatches)

	// filter out matches that are completely contained in others or too short
	rlog.Debugw("culling fragments", "matched", len(featureMatches), "extended", len(extendedMatches))

	// remove extended matches fully enclosed by others
	extendedMatches = cull(extendedMatches, len(feats), 1, 4)

	// create a subject file from the matches' source fragments
	subjectDB, frags := subjectDatabase(extendedMatches, dbs)
	defer os.Remove(subjectDB)

	// re-BLAST the features against the new subject database
	featureMatches = reblastFeatures(identity, feats, conf, subjectDB, frags)

	// merge matches into one another if they can combine to cover a range
	extendedMatches = extendMatches(feats, featureMatches)

	// remove extended matches fully enclosed by others
	extendedMatches = cull(extendedMatches, len(feats), 1, 4)

	rlog.Debugw("culled matches", "remaining", len(extendedMatches))

	// get the full plasmid length as if just synthesizing each feature next to one another
	var targetBuilder strings.Builder
	featureToStart := make(map[int]int) // map from feature index to start index on fullSynthSeq
	for i, feat := range feats {
		featureToStart[i] = targetBuilder.Len()
		targetBuilder.WriteString(feat[1])
	}
	target := targetBuilder.String()

	// get the matches back out of the databases (the full parts)
	frags = []*Frag{}
	seenMatches := make(map[string]bool)
	for _, m := range extendedMatches {
		if _, seen := seenMatches[m.uniqueID]; seen {
			continue
		}
		seenMatches[m.uniqueID] = true

		frag, err := queryDatabases(m.entry, dbs)
		if err != nil {
			rlog.Fatal(err)
		}

		frag.ID = m.entry
		frag.uniqueID = m.uniqueID
		if m.subjectEnd < m.subjectStart {
			continue
		}

		frag.Seq = (frag.Seq + frag.Seq + frag.Seq)[m.subjectStart : m.subjectEnd+1]
		if !m.forward {
			frag.Seq = reverseComplement(frag.Seq)
		}
		frag.conf = conf

		frag.featureStart = m.queryStart
		frag.featureEnd = m.queryEnd

		frag.start = featureToStart[m.queryStart]
		frag.end = featureToStart[(m.queryEnd+1)%len(feats)]

		if frag.end <= frag.start {
			frag.end += len(target) // wrap across the zero index
		}

		if m.queryStart >= len(feats) {
			frag.start += len(target)
			frag.end += len(target)
		}

		frags = append(frags, frag)
	}

	// traverse the fragments, accumulate assemblies that span all the features
	assemblies := createAssemblies(frags, target, len(feats), true, conf)

	// sort assemblies
	sort.Slice(assemblies, func(i, j int) bool {
		return assemblies[i].isBetterThan(assemblies[j])
	})

	var selectedAssemblies []assembly
	if keepNSolutions > 0 {
		if keepNSolutions < len(assemblies) {
			selectedAssemblies = assemblies[:keepNSolutions]
		} else {
			selectedAssemblies = assemblies
		}
	} else {
		// only keep the best solution
		selectedAssemblies = assemblies[:1]
	}

	// fill each assembly and accumulate the pareto optimal solutions
	filledAssemblies := fillAssemblies(target, selectedAssemblies, 0, conf)

	// update the target to the first filled assembly
	if len(filledAssemblies) > 0 {
		target = annealFragments(conf.FragmentsMinHomology, conf.FragmentsMaxHomology, filledAssemblies[0].frags)
	}
	// final sort after filling the assemblies
	sort.Slice(filledAssemblies, func(i, j int) bool {
		return filledAssemblies[i].isBetterThan(*filledAssemblies[j])
	})
	finalSolutions := make([][]*Frag, len(filledAssemblies))
	for i := range finalSolutions {
		finalSolutions[i] = filledAssemblies[i].frags
	}

	return target, finalSolutions
}

// extendMatches groups and extends matches against the subject sequence
func extendMatches(feats [][]string, featureMatches map[string][]featureMatch) (extendedMatches []match) {
	for _, matches := range featureMatches {
		sort.Slice(matches, func(i, j int) bool {
			return matches[i].match.subjectStart < matches[j].match.subjectStart
		})

		// expand the ranges of the matches range based on their continuous feature stretches
		m := matches[0].match
		firstOfStretch := matches[0].featureIndex
		lastOfStretch := matches[0].featureIndex
		numFeatures := 1 // number of features in this extended match

		// features doubled on self to account for feature runs that cross the zero index
		for i, featureMatch := range append(matches, matches...) {
			first := i == 0
			last := i == len(matches)*2-1

			if first {
				continue // still on the first match
			}

			if featureMatch.featureIndex == (lastOfStretch+1)%len(feats) && !last {
				lastOfStretch = featureMatch.featureIndex
				numFeatures++
				if numFeatures < len(feats) {
					continue // continue to extend the feature match stretch
				}
			}

			// cannot extend this stretch, create a new fragment
			m.queryStart = firstOfStretch
			m.queryEnd = lastOfStretch
			m.subjectEnd = matches[(i-1)%len(matches)].match.subjectEnd

			extendedMatches = append(extendedMatches, m)

			// start on the next stretch
			m = featureMatch.match
			firstOfStretch = featureMatch.featureIndex
			lastOfStretch = featureMatch.featureIndex
			numFeatures = 1
		}
	}

	return extendedMatches
}

// create a subject database to query specifically for all
// features. Needed because the first BLAST may not return
// all feature matches on each fragment
func subjectDatabase(extendedMatches []match, dbs []DB) (filename string, frags []*Frag) {
	subject := ""
	for _, m := range extendedMatches {
		frag, err := queryDatabases(m.entry, dbs)
		if err != nil {
			rlog.Fatal(err)
		}
		subject += fmt.Sprintf(">%s\n%s\n", frag.ID, frag.Seq)
		frags = append(frags, frag)
	}

	in, err := os.CreateTemp("", "feature-subject-*")
	if err != nil {
		rlog.Fatal(err)
	}

	if _, err := in.WriteString(subject); err != nil {
		rlog.Fatal(err)
	}

	return in.Name(), frags
}

// reblastFeatures returns matches between the target features and entries in the databases with those features
func reblastFeatures(
	identity int,
	feats [][]string,
	conf *config.Config,
	subjectDB string,
	frags []*Frag) map[string][]featureMatch {
	featureMatches := make(map[string][]featureMatch) // a map from from each entry (by id) to its list of matched features
	for i, target := range feats {
		targetFeature := target[1]
		matches, err := blastAgainst(target[0], targetFeature, subjectDB, identity)
		if err != nil {
			rlog.Fatal(err)
		}

		for _, m := range matches {
			// needs to be at least identity % as long as the queried feature
			mLen := float64(m.subjectEnd - m.subjectStart)
			if mLen/float64(len(targetFeature)) < float64(identity)/100.0 {
				continue
			}

			m.queryStart = i
			m.queryEnd = i
			m.uniqueID = m.entry + strconv.Itoa(m.subjectStart)

			if _, exists := featureMatches[m.entry]; !exists {
				featureMatches[m.entry] = []featureMatch{{featureIndex: i, match: m}}
			} else {
				featureMatches[m.entry] = append(featureMatches[m.entry], featureMatch{featureIndex: i, match: m})
			}
		}
	}

	for i, target := range feats {
		targetFeature := target[1]
		for _, frag := range frags {
			fragSeq := strings.ToUpper(frag.Seq + frag.Seq)
			featureIndex := strings.Index(fragSeq, targetFeature)

			if featureIndex > 0 {
				manualMatch := match{
					entry:        frag.ID,
					uniqueID:     fmt.Sprintf("%s%d", frag.ID, featureIndex),
					querySeq:     targetFeature,
					queryStart:   i,
					queryEnd:     i,
					seq:          targetFeature,
					subjectStart: featureIndex,
					subjectEnd:   featureIndex + len(targetFeature),
					db:           frag.db,
					title:        target[0],
					circular:     frag.fragType == circular,
					forward:      true,
				}

				alreadySeen := false
				for _, m := range featureMatches[frag.ID] {
					if m.featureIndex == i && m.match.subjectStart == manualMatch.subjectStart {
						alreadySeen = true
						break
					}
				}
				if !alreadySeen {
					featureMatches[frag.ID] = append(featureMatches[frag.ID], featureMatch{featureIndex: i, match: manualMatch})
				}
			}
		}
	}

	return featureMatches
}

// NewFeatureDB returns a new copy of the features db
func NewFeatureDB() *kv {
	return newKV(config.FeatureDB)
}

// ListFeatures returns features that are similar in name to the feature name requested.
// if multiple feature names include the feature name, they are all returned.
// otherwise a list of feature names are returned (those beneath a levenshtein distance cutoff)
func ListFeatures(featureName string) {
	f := NewFeatureDB()

	if featureName == "" {
		// no feature name passed, log all of them
		featNames := []string{}
		for feat := range f.contents {
			featNames = append(featNames, feat)
		}

		sort.Slice(
			featNames,
			func(i, j int) bool {
				return strings.ToLower(featNames[i]) < strings.ToLower(featNames[j])
			},
		)

		// print all their names to the console and the first few bp
		w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, '-', tabwriter.TabIndent)
		for _, feat := range featNames {
			seq := f.contents[feat]
			if len(seq) > 20 {
				seq = seq[:20] + "..."
			}
			fmt.Fprintf(w, "%s\t%s\n", feat, seq)
		}

		w.Flush()

		return
	}

	ldCutoff := len(featureName) / 3
	if 1 > ldCutoff {
		ldCutoff = 1
	}
	containing := []string{}
	lowDistance := []string{}

	for fName, fSeq := range f.contents {
		if strings.Contains(fName, featureName) {
			containing = append(containing, fName+"\t"+fSeq)
		} else if len(fName) > ldCutoff && ld(featureName, fName, true) <= ldCutoff {
			lowDistance = append(lowDistance, fName+"\t"+fSeq)
		}
	}

	// check for an exact match
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)
	matchedFeature, exactMatch := f.contents[featureName]
	if exactMatch && len(containing) < 2 {
		fmt.Fprintf(w, featureName+"\t"+matchedFeature)
		if _, err := w.Write([]byte("\n")); err != nil {
			rlog.Fatal(err)
		}
		w.Flush()
		return
	}

	// from https://golang.org/pkg/text/tabwriter/
	if len(containing) < 3 {
		lowDistance = append(lowDistance, containing...)
		containing = []string{} // clear
	}
	if len(containing) > 0 {
		fmt.Fprint(w, strings.Join(containing, "\n"))
	} else if len(lowDistance) > 0 {
		fmt.Fprint(w, strings.Join(lowDistance, "\n"))
	} else {
		if _, err := fmt.Fprintf(w, "failed to find any features for %s", featureName); err != nil {
			rlog.Fatal(err)
		}
	}
	w.Flush()
}

// AddFeatures - add the feature's seq in the database (or create if it isn't in the feature db)
func AddFeatures(name, seq string) {
	f := NewFeatureDB()

	f.contents[name] = seq
	if err := f.save(); err != nil {
		rlog.Fatal(err)
	}
}

// DeleteFeature - delete the feature from the database
func DeleteFeature(name string) {
	f := NewFeatureDB()

	if _, contained := f.contents[name]; !contained {
		fmt.Printf("failed to find %s in the features database\n", name)
	}

	delete(f.contents, name)
	if err := f.save(); err != nil {
		rlog.Fatal(err)
	}
}

// ld compares two strings and returns the levenshtein distance between them.
// This was copied verbatim from https://github.com/spf13/cobra
func ld(s, t string, ignoreCase bool) int {
	if ignoreCase {
		s = strings.ToUpper(s)
		t = strings.ToUpper(t)
	}
	d := make([][]int, len(s)+1)
	for i := range d {
		d[i] = make([]int, len(t)+1)
	}
	for i := range d {
		d[i][0] = i
	}
	for j := range d[0] {
		d[0][j] = j
	}
	for j := 1; j <= len(t); j++ {
		for i := 1; i <= len(s); i++ {
			if s[i-1] == t[j-1] {
				d[i][j] = d[i-1][j-1]
			} else {
				min := d[i-1][j]
				if d[i][j-1] < min {
					min = d[i][j-1]
				}
				if d[i-1][j-1] < min {
					min = d[i-1][j-1]
				}
				d[i][j] = min + 1
			}
		}
	}
	return d[len(s)][len(t)]
}
