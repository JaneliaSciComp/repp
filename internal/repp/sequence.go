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

// SequenceList is for BLAST'ing a sequence against the dbs and finding matches
func SequenceList(
	seq string,
	filters []string,
	identity int,
	leftMargin int,
	dbNames []string) {

	dbs, err := getRegisteredDBs(dbNames)
	if err != nil {
		rlog.Fatal(err)
	}

	matches, err := blast("find_cmd", seq, true, leftMargin, dbs, filters, identity, false)
	if err != nil {
		rlog.Fatal(err)
	}

	if len(matches) == 0 {
		rlog.Fatal("no matches found")
	}

	// sort so the largest matches are first
	sort.Slice(matches, func(i, j int) bool {
		return (matches[i].subjectEnd - matches[i].subjectStart) > (matches[j].queryEnd - matches[j].queryStart)
	})

	// to avoid logging the same matches multiple times
	key := func(m match) string {
		return m.entry + strconv.Itoa(m.subjectStart) + strconv.Itoa(m.subjectEnd)
	}

	seenIds := make(map[string]bool)
	writer := tabwriter.NewWriter(os.Stdout, 0, 4, 3, ' ', 0)
	fmt.Fprintf(writer, "entry\tqstart\tqend\tsstart\tsend\tdatabase\t\n")
	for _, m := range matches {
		if _, seen := seenIds[key(m)]; seen {
			continue
		}

		if m.subjectEnd-m.subjectStart < 20 {
			continue
		}

		fmt.Fprintf(writer, "%s\t%d\t%d\t%d\t%d\t%s\n", m.entry, m.queryStart, m.queryEnd, m.subjectStart, m.subjectEnd, m.db.Name)
		seenIds[key(m)] = true
	}
	writer.Flush()
}

// Sequence is for running an end to end plasmid design using a target sequence.
func Sequence(assemblyParams AssemblyParams, maxSolutions int, conf *config.Config) (solutions [][]*Frag) {
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
	// build up the assemblies that make the sequence
	target, solutions, err := sequence(
		assemblyParams.GetIn(),
		assemblyParams.GetFilters(),
		assemblyParams.GetIdentity(),
		assemblyParams.GetLeftMargin(),
		backboneFrag,
		dbs,
		maxSolutions,
		conf)
	if err != nil {
		rlog.Fatal(err)
	}

	primersDB := readOligos(assemblyParams.GetPrimersDBLocations(), primerIDPrefix, false)
	synthFragsDB := readOligos(assemblyParams.GetSynthFragsDBLocations(), synthFragIDPrefix, true)

	// write the results to a file
	elapsed := time.Since(start)
	_, err = writeResult(
		assemblyParams.GetOut(),
		assemblyParams.GetOutputFormat(),
		target.ID,
		target.Seq,
		solutions,
		primersDB,
		synthFragsDB,
		backboneMeta,
		elapsed.Seconds(),
		conf,
	)
	if err != nil {
		rlog.Fatal(err)
	}

	rlog.Debugw("execution time", "execution", elapsed)

	return solutions
}

// sequence builds a plasmid cost optimization
//
// The goal is to find an "optimal" assembly sequence with:
//  1. the fewest fragments
//  2. the lowest overall assembly cost ($)
//
// and, secondarily:
//  3. no duplicate end regions between Gibson fragments
//  4. no hairpins in the junctions
//  5. no off-target binding sites in the parent plasmids
//  6. low primer3 penalty scores
//
// First build up assemblies, creating all possible assemblies that are
// beneath the upper-bound limit on the number of fragments fully covering
// the target sequence
//
// Then find the pareto optimal solutions that minimize either the cost
// of the assembly or the number of fragments (relative to the other
// assembly plans)
//
// Then, for each pareto optimal solution, traverse the assembly and
// "fill-in" the nodes. Create primers on the Frag if it's a PCR Frag
// or create a sequence to be synthesized if it's a synthetic fragment.
// Error out and repeat the build stage if a Frag fails to be filled
func sequence(
	input string,
	filters []string,
	identity int,
	leftMargin int,
	backboneFrag *Frag,
	dbs []DB,
	keepNSolutions int,
	conf *config.Config) (target *Frag, solutions [][]*Frag, err error) {

	// read the target sequence (the first in the slice is used)
	fragments, err := read(input, false, false)
	if err != nil {
		return &Frag{}, nil, fmt.Errorf("failed to read target sequence from %s: %v", input, err)
	}

	if len(fragments) > 1 {
		rlog.Warnf(
			"warning: %d fragments were in %s. Only targeting the sequence of the first: %s\n",
			len(fragments),
			input,
			fragments[0].ID,
		)
	}

	target = fragments[0]
	targetSeqLen := len(target.Seq)
	rlog.Debugw("building plasmid", "targetID", target.ID, "targetLen", targetSeqLen)

	var bbFragInsert *Frag
	if backboneFrag.ID != "" {
		bbSeqLen := len(backboneFrag.Seq)
		inputSeq := strings.ToUpper(target.Seq + target.Seq)
		bbSubSeqIndex := strings.Index(inputSeq, backboneFrag.Seq)
		if bbSubSeqIndex != -1 {
			// add the backbone to the sequence of the target plasmid
			bbFragInsert = &Frag{
				ID:         backboneFrag.ID,
				Seq:        backboneFrag.Seq,
				uniqueID:   "backbone" + strconv.Itoa(bbSubSeqIndex),
				start:      bbSubSeqIndex,
				end:        bbSubSeqIndex + bbSeqLen,
				fragType:   pcr,
				matchRatio: 1,
				conf:       conf,
			}
		} else {
			revBBSeq := reverseComplement(backboneFrag.Seq)
			revBBSubSeqIndex := strings.Index(inputSeq, revBBSeq)
			if revBBSubSeqIndex != -1 {
				// add the reverse backbone to the sequence of the target plasmid
				bbFragInsert = &Frag{
					ID:         backboneFrag.ID,
					Seq:        revBBSeq,
					uniqueID:   "backbone" + strconv.Itoa(revBBSubSeqIndex),
					start:      revBBSubSeqIndex,
					end:        revBBSubSeqIndex + bbSeqLen,
					fragType:   pcr,
					matchRatio: 1,
					conf:       conf,
				}
			} else {
				// add the backbone to the sequence of the target plasmid
				bbFragInsert = &Frag{
					ID:         backboneFrag.ID,
					Seq:        backboneFrag.Seq,
					uniqueID:   "backbone" + strconv.Itoa(targetSeqLen),
					start:      targetSeqLen,
					end:        targetSeqLen + bbSeqLen,
					fragType:   pcr,
					matchRatio: 1,
					conf:       conf,
				}
				target.Seq += backboneFrag.Seq
			}
		}
	} else {
		bbFragInsert = nil
	}

	// get all the matches against the target plasmid
	matches, err := blast(
		target.ID,
		target.Seq,
		true,
		leftMargin,
		dbs,
		filters,
		identity,
		conf.UseUngappedAlignment,
	)
	if err != nil {
		dbMessage := strings.Join(dbNames(dbs), ", ")
		return &Frag{}, nil, fmt.Errorf("failed to blast %s against the dbs %s: %v", target.ID, dbMessage, err)
	}

	// keep only "proper" arcs (non-self-contained)
	matches = cull(matches, conf.PcrMinFragLength, 1)
	rlog.Debugw("culled matches", "remaining", len(matches)/2)

	// map fragment Matches to nodes
	frags := newFrags(matches, conf)

	if bbFragInsert != nil {
		copiedBB := bbFragInsert.copy()
		copiedBB.start += len(target.Seq)
		copiedBB.end += len(target.Seq)
		// add the backbone in as fragment (copy twice across zero index)
		frags = append(frags, bbFragInsert, copiedBB)
		sort.Slice(frags, func(i, j int) bool {
			return frags[i].start < frags[j].start
		})
	}

	// build up a slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target plasmid
	assemblies := createAssemblies(frags, target.Seq, len(target.Seq), false, conf)

	rlog.Debugf("Sort %d found assemblies\n", len(assemblies))
	// sort assemblies
	sort.Slice(assemblies, func(i, j int) bool {
		return assemblies[i].isBetterThan(assemblies[j])
	})
	if isVerboseLogging() {
		for i, a := range assemblies {
			rlog.Debugf("Prelim solution %d: %v", i+1, a)
		}
	}
	var maxSolutions int
	if keepNSolutions > 0 {
		if keepNSolutions < len(assemblies) {
			maxSolutions = keepNSolutions
		} else {
			maxSolutions = len(assemblies)
		}
	} else {
		// only keep the best solution
		maxSolutions = 1
	}
	maxInspectedSolutions := maxSolutions + int(0.2*float32(len(assemblies)))

	var filledAssemblies []*assembly

	rlog.Infof("Start filling PCR primers for %d assemblies out of %d\n", maxSolutions, len(assemblies))
	// try to fill as many solutions as requested (if there are enough assemblies)
	// so if not all solutions could be filled try other assemblies
	for searchSolutionFromIndex := 0; searchSolutionFromIndex < len(assemblies); searchSolutionFromIndex += maxInspectedSolutions {
		var selectedAssemblies []assembly
		var lastInspectedIndex = searchSolutionFromIndex + maxInspectedSolutions - len(filledAssemblies)
		if lastInspectedIndex < len(assemblies) {
			rlog.Infof("Inspecting and filling assemblies from %d to %d", searchSolutionFromIndex, lastInspectedIndex)
			selectedAssemblies = assemblies[searchSolutionFromIndex:lastInspectedIndex]
		} else {
			rlog.Infof("Inspecting and filling assemblies from %d to the end", searchSolutionFromIndex)
			selectedAssemblies = assemblies[searchSolutionFromIndex:]
		}
		// fill in only top best assemblies
		solutions := fillAssemblies(target.Seq, selectedAssemblies, searchSolutionFromIndex, conf)
		filledAssemblies = append(filledAssemblies, solutions...)
		if len(filledAssemblies) >= maxSolutions {
			break
		} else {
			rlog.Infof("Filled %d solutions out of the first %d assemblies\n",
				len(filledAssemblies),
				searchSolutionFromIndex+len(selectedAssemblies))
			if searchSolutionFromIndex+len(selectedAssemblies) < len(assemblies) {
				rlog.Infof("Try to fill remaining %d solutions out of %d found assemblies\n",
					maxSolutions-len(filledAssemblies),
					len(assemblies)-searchSolutionFromIndex-len(selectedAssemblies))
			}
		}
	}
	// final sort after filling the assemblies
	// but this time sort by the number of fragments
	sort.Slice(filledAssemblies, func(i, j int) bool {
		return filledAssemblies[i].len() < filledAssemblies[j].len()
	})
	rlog.Infof("Finished filling %d assemblies", len(filledAssemblies))
	var nfinalSolutions int
	if len(filledAssemblies) < maxSolutions {
		nfinalSolutions = len(filledAssemblies)
	} else {
		nfinalSolutions = maxSolutions
	}
	finalSolutions := make([][]*Frag, nfinalSolutions)
	for i := range finalSolutions {
		finalSolutions[i] = filledAssemblies[i].frags
	}
	return target, finalSolutions, nil
}
