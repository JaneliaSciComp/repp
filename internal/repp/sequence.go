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
	dbNames []string) {

	dbs, err := getRegisteredDBs(dbNames)
	if err != nil {
		rlog.Fatal(err)
	}

	matches, err := blast("find_cmd", seq, true, dbs, filters, identity)
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
func Sequence(assemblyParams AssemblyParams, conf *config.Config) (solutions [][]*Frag) {
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
	insert, target, solutions, err := sequence(
		assemblyParams.GetIn(),
		assemblyParams.GetFilters(),
		assemblyParams.GetIdentity(),
		backboneFrag,
		dbs,
		conf)
	if err != nil {
		rlog.Fatal(err)
	}

	// write the results to a file
	elapsed := time.Since(start)
	_, err = writeJSON(
		assemblyParams.GetOut(),
		target.ID,
		target.Seq,
		solutions,
		len(insert.Seq),
		elapsed.Seconds(),
		backboneMeta,
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
	backboneFrag *Frag,
	dbs []DB,
	conf *config.Config) (insert, target *Frag, solutions [][]*Frag, err error) {

	// read the target sequence (the first in the slice is used)
	fragments, err := read(input, false)
	if err != nil {
		return &Frag{}, &Frag{}, nil, fmt.Errorf("failed to read target sequence from %s: %v", input, err)
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
	rlog.Debugw("building plasmid", "targetID", target.ID)

	// if a backbone was specified, add it to the sequence of the target frag
	insert = target.copy() // store a copy for logging later
	if backboneFrag.ID != "" {
		target.Seq += backboneFrag.Seq
	}

	// get all the matches against the target plasmid
	matches, err := blast(target.ID, target.Seq, true, dbs, filters, identity)
	if err != nil {
		dbMessage := strings.Join(dbNames(dbs), ", ")
		return &Frag{}, &Frag{}, nil, fmt.Errorf("failed to blast %s against the dbs %s: %v", target.ID, dbMessage, err)
	}

	// keep only "proper" arcs (non-self-contained)
	matches = cull(matches, len(target.Seq), conf.PcrMinLength, 1)
	rlog.Debugw("culled matches", "remaining", len(matches)/2)

	// map fragment Matches to nodes
	frags := newFrags(matches, conf)

	if backboneFrag.ID != "" {
		// add the backbone in as fragment (copy twice across zero index)
		backboneFrag.conf = conf
		backboneFrag.start = len(insert.Seq)
		backboneFrag.end = backboneFrag.start + len(backboneFrag.Seq) - 1
		backboneFrag.uniqueID = "backbone" + strconv.Itoa(backboneFrag.start)
		frags = append(frags, backboneFrag)

		copiedBB := backboneFrag.copy()
		copiedBB.start += len(target.Seq)
		copiedBB.end += len(target.Seq)
		copiedBB.uniqueID = backboneFrag.uniqueID
		frags = append(frags, copiedBB)

		sort.Slice(frags, func(i, j int) bool {
			return frags[i].start < frags[j].start
		})
	}

	// build up a slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target plasmid
	assemblies := createAssemblies(frags, target.Seq, len(target.Seq), false, conf)

	// build up a map from fragment count to a sorted list of assemblies with that number
	assemblyCounts, countToAssemblies := groupAssembliesByCount(assemblies)

	// fill in pareto optimal assembly solutions
	solutions = fillAssemblies(target.Seq, assemblyCounts, countToAssemblies, conf)

	return insert, target, solutions, nil
}
