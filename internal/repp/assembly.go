package repp

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target plasmid.
type assembly struct {
	// frags, ordered by distance from the "end" of the plasmid
	frags []*Frag

	// estimated cost of making this assembly
	cost float64

	// cost adjusted based on synthetic fragments
	adjustedCost float64

	// total number of synthetic nodes that will be needed to make this
	synths int
}

// createNewAssembly Frag to the end of an assembly. Return a new assembly and whether it circularized
func createNewAssembly(existingAssembly assembly, f *Frag, maxCount, targetLength int,
	features bool) (assembly /*created*/, bool /*circularized*/, bool) {

	first := existingAssembly.firstFrag()
	last := existingAssembly.lastFrag()

	var existingAssemblyStart int
	var existingAssemblyEnd int
	var start int
	var end int

	if features {
		existingAssemblyStart = first.featureStart
		existingAssemblyEnd = last.featureEnd
		start = f.featureStart
		end = f.featureEnd
	} else {
		existingAssemblyStart = first.start
		existingAssemblyEnd = last.end
		start = f.start
		end = f.end
	}

	// check if we could complete an assembly with this new Frag
	circularized := end >= existingAssemblyStart+targetLength-1

	// check if this is the first fragment annealing to itself
	selfAnnealing := f.uniqueID == first.uniqueID

	// calc the number of synthesis fragments needed to get to this next Frag
	synths := last.synthDist(f)
	if features && start > existingAssemblyEnd {
		synths = start - existingAssemblyEnd - 1
	}

	newCount := existingAssembly.len() + synths
	if !selfAnnealing {
		newCount++
	}

	assemblyEnd := existingAssemblyEnd
	if newCount > maxCount || (end-assemblyEnd < f.conf.PcrMinLength && !features) {
		return assembly{}, false, false
	}

	// calc the estimated dollar cost of getting to the next Frag
	annealCost, adjustedCost := last.costTo(f)
	if selfAnnealing && synths == 0 {
		annealCost = 0   // does not cost extra to anneal to the first fragment
		adjustedCost = 0 // there are no synth so it is safe to set it to 0
	}

	// check whether the Frag is already contained in the assembly
	// if so, the cost of procurement is not incurred twice
	fragContained := false
	for _, included := range existingAssembly.frags {
		if included.ID == f.ID && included.fragType == f.fragType {
			fragContained = true
			break
		}
	}

	if fragContained {
		// don't double count the cost of procuring this Frag to the total assembly cost
		fragCost, adjustedFragCost := f.cost(false)
		annealCost += fragCost
		adjustedCost += adjustedFragCost
	} else {
		fragCost, adjustedFragCost := f.cost(true)
		annealCost += fragCost
		adjustedCost += adjustedFragCost
	}

	// copy over all the fragments, need to avoid referencing same frags
	newFrags := []*Frag{}
	for _, frag := range existingAssembly.frags {
		newFrags = append(newFrags, frag.copy())
	}
	if !selfAnnealing {
		newFrags = append(newFrags, f.copy())
	}

	return assembly{
		frags:        newFrags,
		cost:         existingAssembly.cost + annealCost,
		adjustedCost: existingAssembly.adjustedCost + adjustedCost,
		synths:       existingAssembly.synths + synths,
	}, true, circularized
}

func (a assembly) firstFrag() *Frag {
	return a.frags[0]
}

func (a assembly) lastFrag() *Frag {
	return a.frags[len(a.frags)-1]
}

// len returns len(assembly.nodes) + the synthesis fragment count.
func (a assembly) len() int {
	return len(a.frags) + a.synths
}

// fill traverses frags in an assembly and adds primers or makes synthetic fragments where necessary.
// It can fail. For example, a PCR Frag may have off-targets in the parent plasmid.
func (a assembly) fill(target string, conf *config.Config) ([]*Frag, error) {
	// check for and error out if there are duplicate ends between fragments,
	// ie unintended junctions between fragments that shouldn't be annealing
	if hasDuplicate, left, right, dupSeq := duplicates(a.frags, conf.FragmentsMinHomology, conf.FragmentsMaxHomology); hasDuplicate {
		return nil, fmt.Errorf("duplicate junction between %s and %s: %s", left, right, dupSeq)
	}

	// edge case where a single Frag fills the whole target plasmid. Return just a single
	// "fragment" (of circular type... it is misnomer) that matches the target sequence 100%
	if a.len() == 1 && len(a.frags[0].Seq) >= len(target) {
		f := a.frags[0]

		return []*Frag{
			{
				ID:         f.ID,
				Seq:        strings.ToUpper(f.Seq)[0:len(target)], // it may be longer
				fragType:   circular,
				matchRatio: f.matchRatio,
				conf:       conf,
			},
		}, nil
	}

	// copy all the fragments. needed because ranges are mutated in assembly.fill,
	// so distance to neightbor estimates become invalid after a neighbor is mutated
	var origFrags []*Frag
	for _, f := range a.frags {
		origFrags = append(origFrags, f.copy())
	}

	pcrFrags := []*Frag{}

	// fill in primers. let each Frag create primers for itself that
	// will span it to the last and next fragments (if reachable)
	for i, f := range a.frags {
		// try and make primers for the fragment (need prev and next nodes)
		prev := prevFragment(origFrags, i, target, conf)
		next := nextFragment(origFrags, i, target, conf)

		needsPCR := f.fragType == circular ||
			f.fragType == pcr ||
			!prev.overlapsViaHomology(f) && prev.couldOverlapViaPCR(f) ||
			!f.overlapsViaHomology(next) && f.couldOverlapViaPCR(next)

		// if the Frag has a full target from upload or
		if needsPCR {
			// create primers for the Frag and add them to the Frag if it needs them
			// to anneal to the adjacent fragments
			if err := f.setPrimers(prev, next, target, conf); err != nil || len(f.Primers) < 2 {
				return nil, fmt.Errorf("failed to pcr %s: %v", f.ID, err)
			}
			f.fragType = pcr // is now a pcr type
		}

		// accumulate the prepared fragment
		pcrFrags = append(pcrFrags, f)
	}

	// second loop to fill in gaps between fragments that need to be filled via synthesis
	pcrAndSynthFrags := []*Frag{}
	for i, f := range pcrFrags {
		if f.Seq != "" {
			pcrAndSynthFrags = append(pcrAndSynthFrags, f)
		}

		// add synthesized fragments between this Frag and the next (if necessary)
		next := nextFragment(pcrFrags, i, target, conf)
		if synthedFrags := f.synthTo(next, target); synthedFrags != nil {
			pcrAndSynthFrags = append(pcrAndSynthFrags, synthedFrags...)
		}
	}
	// validate that fragments will anneal to one another
	if err := validateJunctions(pcrAndSynthFrags, conf); err != nil {
		return nil, err
	}

	return pcrAndSynthFrags, nil
}

// nextFragment returns the fragment that's one beyond the one passed.
// The fragments are considered to be part of a "circular" sequence
// simulated by concatenating the sequence to itself
// the next fragment after the last from the list is based on the
// first fragment from the list by adding target sequence length to its start and end
func nextFragment(frags []*Frag, i int, target string, conf *config.Config) *Frag {
	if i < len(frags)-1 {
		return frags[i+1]
	}

	// mock up a next fragment that's to the right of this terminal Frag
	return &Frag{
		start: frags[0].start + len(target),
		end:   frags[0].end + len(target),
		conf:  conf,
	}
}

// prevFragment returns the fragment that's one before the current one.
// The fragments are considered to be part of a "circular" sequence
// simulated by concatenating the sequence to itself
// the prev fragment of the first from the list is based on the
// last fragment from the list by subtracting the length of the target sequence
// from its start and end
func prevFragment(frags []*Frag, i int, target string, conf *config.Config) *Frag {
	if i > 0 {
		return frags[i-1]
	}

	// mock up a next fragment that's to the right of this terminal Frag
	return &Frag{
		start: frags[len(frags)-1].start - len(target),
		end:   frags[len(frags)-1].end - len(target),
		conf:  conf,
	}
}

// duplicates runs through all the nodes in an assembly and checks whether any of
// them have unintended homology, or "duplicate homology".
func duplicates(frags []*Frag, min, max int) (isDup bool, first, second, dup string) {
	c := len(frags) // Frag count
	for i, f := range frags {
		// check to make sure the fragment doesn't anneal to itself
		if c > 1 {
			if selfJ := f.selfJunction(min, max); selfJ != "" && len(selfJ) < len(f.Seq) {
				return true, f.ID, f.ID, selfJ
			}
		}

		for j := 2; j < c; j++ { // skip next Frag, i+1 is supposed to anneal to i
			junc := f.junction(frags[(j+i)%c], min, max)
			if junc != "" {
				return true, f.ID, frags[(j+i)%c].ID, junc
			}
		}
	}

	return false, "", "", ""
}

// createAssemblies builds up circular assemblies (unfilled lists of fragments that should be combinable)
//
// It is created by traversing a DAG in forward order:
//
// foreach fragment (sorted in increasing start index order):
//
//	  foreach otherFragment that fragment overlaps with + reachSynthCount more:
//		   foreach assembly on fragment:
//	      add otherFragment to the assembly to create a new assembly, store on otherFragment
func createAssemblies(frags []*Frag, target string, targetLength int, features bool, conf *config.Config) []assembly {
	// sort by start index again
	sort.Slice(frags, func(i, j int) bool {
		return frags[i].start < frags[j].start
	})

	// indexedAssemblies[i] holds all assemblies that are extended starting with the i-th fragment
	var indexedAssemblies = make([][]assembly, len(frags))

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target plasmid... 100% match
		// it is the target plasmid. just return that as the assembly
		if len(f.Seq) >= targetLength && !features {
			return []assembly{
				{
					frags:  []*Frag{f.copy()},
					synths: 0,
				},
			}
		}

		// create a starting assembly for each fragment containing just it
		cost, adjustedCost := f.costTo(f)
		indexedAssemblies[i] = []assembly{
			{
				frags:        []*Frag{f.copy()}, // just self
				cost:         cost,              // just PCR,
				adjustedCost: adjustedCost,
				synths:       0, // no synthetic frags at start
			},
		}
	}

	finalAssemblies := []assembly{}

	for i, f := range frags { // for every Frag in the list of increasing start index frags
		for _, j := range f.reach(frags, i, features) { // for every overlapping fragment + reach more
			for _, a := range indexedAssemblies[i] { // for every assembly on the reaching fragment
				newAssembly, created, circularized := createNewAssembly(a, frags[j], conf.FragmentsMaxCount, targetLength, features)

				if !created { // if a new assembly wasn't created, move on
					continue
				}

				if circularized { // we've circularized a plasmid, it's ready for filling
					finalAssemblies = append(finalAssemblies, newAssembly)
				} else {
					// the new fragment was created by adding the j-th fragment
					// so when it's processed it will be extended started with j-th frag
					// this works because j > i so indexedAssemblies[j] is still in the queue
					indexedAssemblies[j] = append(indexedAssemblies[j], newAssembly)
				}
			}
		}
	}

	// create a fully synthetic plasmid from just synthetic fragments
	// in case all other plasmid designs fail
	mockStart := &Frag{start: conf.FragmentsMinHomology, end: conf.FragmentsMinHomology, conf: conf}
	mockEnd := &Frag{start: len(target), end: len(target), conf: conf}
	cost, adjustedCost := mockStart.costTo(mockEnd)
	synths := mockStart.synthTo(mockEnd, target)
	finalAssemblies = append(finalAssemblies, assembly{
		frags:        synths,
		cost:         cost,
		adjustedCost: adjustedCost,
		synths:       len(synths),
	})
	rlog.Debugw("assemblies made", "count", len(finalAssemblies))

	return finalAssemblies
}

// groupAssembliesByCount returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost.
func groupAssembliesByCount(assemblies []assembly) ([]int, map[int][]assembly) {
	countToAssemblies := make(map[int][]assembly)
	for _, a := range assemblies {
		if as, ok := countToAssemblies[a.len()]; ok {
			countToAssemblies[a.len()] = append(as, a)
		} else {
			countToAssemblies[a.len()] = []assembly{a}
		}
	}

	// sort the fragment counts of assemblies and the assemblies within each
	// assembly count, so we're trying the shortest assemblies first, and the cheapest
	// assembly within each fragment count before the others
	var counts []int
	for count := range countToAssemblies {
		counts = append(counts, count)
		sort.Slice(countToAssemblies[count], func(i, j int) bool {
			return countToAssemblies[count][i].cost < countToAssemblies[count][j].cost
		})
	}
	sort.Ints(counts)

	return counts, countToAssemblies
}

// fillAssemblies fills in assemblies and returns the pareto optimal solutions.
func fillAssemblies(target string, counts []int, countToAssemblies map[int][]assembly, conf *config.Config) (solutions [][]*Frag) {
	// append a fully synthetic solution at first, nothing added should cost more than this (single plasmid)
	filled := make(map[int][]*Frag)
	minCostAssembly := math.MaxFloat64

	for _, count := range counts {
		for _, assemblyToFill := range countToAssemblies[count] {
			if assemblyToFill.adjustedCost > minCostAssembly {
				// skip this and the rest with this count, there's another
				// cheaper option with the same number or fewer fragments (estimated)
				break
			}

			filledFragments, err := assemblyToFill.fill(target, conf)
			if err != nil || filledFragments == nil {
				// assemblyToFill.log()
				// fmt.Println("error", err.Error())
				// Log.Fatal(err)
				continue
			}

			_, adjustedAssemblyCost := fragsCost(filledFragments)

			if adjustedAssemblyCost >= minCostAssembly || len(filledFragments) > conf.FragmentsMaxCount {
				continue // wasn't actually cheaper, keep trying
			}
			minCostAssembly = adjustedAssemblyCost // store this as the new cheapest assembly

			// delete all assemblies with more fragments that cost more
			for filledCount, existingFilledFragments := range filled {
				if filledCount < len(filledFragments) {
					continue
				}

				_, existingAdjustedCost := fragsCost(existingFilledFragments)
				if existingAdjustedCost >= adjustedAssemblyCost {
					delete(filled, filledCount)
				}
			}

			// set this is as the new cheapest of this length
			filled[len(filledFragments)] = filledFragments
		}
	}

	for _, frags := range filled {
		solutions = append(solutions, frags) // flatten
	}

	return solutions
}
