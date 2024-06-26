package repp

import (
	"fmt"
	"sort"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
	"golang.org/x/exp/maps"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target plasmid.
type assembly struct {
	// frags, ordered by distance from the "end" of the plasmid
	frags []*Frag

	// self annealed - last and first fragment are identical
	selfAnnealing bool

	// estimated cost of making this assembly
	cost float64

	// cost adjusted based on synthetic fragments
	adjustedCost float64

	// number of pcr frags
	pcrs int
	// total number of synthetic nodes that will be needed to make this
	synths int
}

// String display method for an assembly
func (a assembly) String() string {
	res := ""
	for i, f := range a.frags {
		if i > 0 {
			res = res + " -> " + f.String()
		} else {
			if a.selfAnnealing {
				res = "(+)" + f.String()
			} else {
				res = f.String()
			}
		}
	}
	return fmt.Sprintf("%s (n=%d, nfrags=%d, coverage=%d, c=%f, ac=%f)",
		res, a.len(), len(a.frags), a.coverage(),
		a.cost, a.adjustedCost)
}

// return assembly hash based on fragment IDs
func (a assembly) assemblyHash() string {
	fragIDs := map[string]int8{}
	for _, f := range a.frags {
		fragIDs[f.uniqueID] = 1
	}
	return fmt.Sprintf("%v", fragIDs)
}

// get the first fragment of the assembly
func (a assembly) firstFrag() *Frag {
	return a.frags[0]
}

// get the last fragment of the assembly
func (a assembly) lastFrag() *Frag {
	return a.frags[len(a.frags)-1]
}

// len returns len(assembly.nodes) + the synthesis fragment count.
func (a assembly) len() int {
	return len(a.frags) + a.synths
}

// get the count of covered bps
func (a assembly) coverage() int {
	bps := 0
	for _, f := range a.frags {
		if f.matchRatio > 0 {
			bps += len(f.Seq)
		}
	}
	return bps
}

func (a assembly) isBetterThan(ref assembly) bool {
	if a.len() < ref.len() {
		return true
	} else if a.len() > ref.len() {
		return false
	}
	if a.synths < ref.synths {
		return true
	} else if a.synths > ref.synths {
		return false
	}
	return a.adjustedCost <= ref.adjustedCost
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
				uniqueID:   f.uniqueID,
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
				return nil, err
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
		return pcrAndSynthFrags, err
	}

	return pcrAndSynthFrags, nil
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
	rlog.Debugf("Fragments selected to create the assembly: %v\n", frags)
	// indexedAssemblies[i] holds all assemblies that are extended
	// from the i-th element of frags
	var indexedAssemblies = make([][]assembly, len(frags))

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target plasmid... 100% match
		// it is the target plasmid. just return that as the assembly
		if len(f.Seq) >= targetLength && !features {
			rlog.Infof("Target completelly covered by a single plasmid assembly")
			return []assembly{
				{
					frags:  []*Frag{f.copy()},
					synths: 0,
					pcrs:   1,
				},
			}
		}
		// create a starting assembly for each fragment containing just it
		cost, adjustedCost := f.cost(true)
		indexedAssemblies[i] = []assembly{
			{
				frags:        []*Frag{f.copy()}, // just self
				cost:         cost,              // just PCR,
				adjustedCost: adjustedCost,
				synths:       0, // no synthetic frags at start
				pcrs:         1,
			},
		}
	}

	finalAssemblies := map[string]assembly{}

	for i, f := range frags { // for every Frag in the list of increasing start index frags
		for _, j := range f.reach(frags, i, features) { // for every overlapping fragment + reach more
			for _, a := range indexedAssemblies[i] { // for every assembly on the reaching fragment
				rlog.Debugf("Trying to extend %v with %v", a, frags[j])
				newAssembly, complete, err := extendAssembly(a, frags[j], conf.FragmentsMaxCount, targetLength, features)
				if err != nil { // if a new assembly wasn't created, move on
					rlog.Debugf("%v could not be extended with %v because %v", a, frags[j], err)
					continue
				}

				if complete { // we've circularized a plasmid, it's ready for filling
					newAssemblyID := newAssembly.assemblyHash()
					if _, exists := finalAssemblies[newAssemblyID]; !exists {
						rlog.Debugf("Adding final assembly: %v", newAssembly)
						finalAssemblies[newAssemblyID] = newAssembly
					} else {
						rlog.Debugf("Discard %v - was already found", newAssembly)
					}
				} else {
					// the new fragment was created by adding the j-th fragment
					// so when it's processed it will be extended started with j-th frag
					// this works because j > i so indexedAssemblies[j] is still in the queue

					// before considering it check that it has not already reached the allowed number of fragments
					if newAssembly.len() < conf.FragmentsMaxCount {
						indexedAssemblies[j] = append(indexedAssemblies[j], newAssembly)
					} else {
						// if a is already at the max length and it's not complete so do not even attempt to extend this anymore
						rlog.Debugf("Abandon candidate %v because it already reached the max fragments count: %d\n",
							newAssembly, newAssembly.len())
						continue
					}
				}
			}
		}
	}

	// create a fully synthetic plasmid from just synthetic fragments
	// in case all other plasmid designs fail
	mockStart := &Frag{
		uniqueID: "mockStart",
		start:    conf.FragmentsMinHomology,
		end:      conf.FragmentsMinHomology,
		conf:     conf,
	}
	mockEnd := &Frag{
		uniqueID: "mockEnd",
		start:    len(target),
		end:      len(target),
		conf:     conf,
	}
	cost, adjustedCost := mockStart.costTo(mockEnd)
	synths := mockStart.synthTo(mockEnd, target)
	mockSynthAssembly := assembly{
		frags:        synths,
		cost:         cost,
		adjustedCost: adjustedCost,
		synths:       len(synths),
		pcrs:         0,
	}
	if _, mockAssemblyFound := finalAssemblies[mockSynthAssembly.assemblyHash()]; mockAssemblyFound {
		rlog.Errorf("Found an assembly similar to the mock synthesized assembly: %v", mockSynthAssembly)
	} else {
		finalAssemblies[mockSynthAssembly.assemblyHash()] = mockSynthAssembly
	}
	rlog.Infof("Found a total of %d assemblies", len(finalAssemblies))

	return maps.Values(finalAssemblies)
}

// extendAssembly - extends currentAssembly by add a new Frag to its end.
// Return the new extended assembly and whether it is complete
func extendAssembly(currentAssembly assembly, f *Frag, maxCount, targetLength int,
	features bool) (assembly, bool, error) {

	first := currentAssembly.firstFrag()
	last := currentAssembly.lastFrag()

	var currentAssemblyStart int
	var currentAssemblyEnd int
	var start int
	var end int

	if features {
		currentAssemblyStart = first.featureStart
		currentAssemblyEnd = last.featureEnd
		start = f.featureStart
		end = f.featureEnd
	} else {
		currentAssemblyStart = first.start
		currentAssemblyEnd = last.end
		start = f.start
		end = f.end
	}

	// check if we could complete an assembly with this new Frag
	complete := end >= currentAssemblyStart+targetLength-1

	// check if this is the first fragment annealing to itself
	selfAnnealing := f.uniqueID == first.uniqueID

	// calc the number of synthesis fragments needed to get to this next Frag
	synths := last.synthDist(f)
	if features && start > currentAssemblyEnd {
		synths = start - currentAssemblyEnd - 1
	}

	newCount := currentAssembly.len() + synths
	if !selfAnnealing {
		newCount++
	}

	assemblyEnd := currentAssemblyEnd
	if newCount > maxCount {
		return assembly{}, false, fmt.Errorf("the resulted assembly has  more fragments than allowed (%d > %d)", newCount, maxCount)
	}
	if end-assemblyEnd < f.conf.PcrMinFragLength && !features {
		return assembly{}, false, fmt.Errorf("overlap with last fragment is too short (%d < %d)", end-assemblyEnd, f.conf.PcrMinFragLength)
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
	for _, included := range currentAssembly.frags {
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
	for _, frag := range currentAssembly.frags {
		newFrags = append(newFrags, frag.copy())
	}
	if !selfAnnealing {
		newFrags = append(newFrags, f.copy())
	}

	return assembly{
		frags:         newFrags,
		selfAnnealing: selfAnnealing,
		cost:          currentAssembly.cost + annealCost,
		adjustedCost:  currentAssembly.adjustedCost + adjustedCost,
		synths:        currentAssembly.synths + synths,
		pcrs:          currentAssembly.pcrs + 1,
	}, complete, nil
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

// fillAssemblies fills in assemblies and returns the pareto optimal solutions.
func fillAssemblies(target string, assemblies []assembly, selectedAssembliesStart int, conf *config.Config) (solutions []*assembly) {
	var filled []*assembly
	for ai, a := range assemblies {
		rlog.Debugf("Try to fill a[%d]: %v\n", selectedAssembliesStart+ai+1, a)
		filledFragments, err := a.fill(target, conf)
		if err != nil || filledFragments == nil || len(filledFragments) == 0 {
			// this error can be pretty verbose so I am only displaying it in debug mode
			rlog.Debugf("Error filling assembly a[%d]: %v because: %v\n",
				selectedAssembliesStart+ai+1, a, err)
		} else {
			assemblyCost := 0.0
			assemblyAdjustedCost := 0.0
			npcrs := 0
			nsynths := 0
			for _, f := range filledFragments {
				if f.fragType == pcr {
					npcrs++
				} else {
					nsynths++
				}
				// assume no procurement cost
				fCost, fAdjustedCost := f.cost(false)
				assemblyCost += fCost
				assemblyAdjustedCost += fAdjustedCost
			}
			filledAssembly := &assembly{
				frags:        filledFragments,
				cost:         assemblyCost,
				adjustedCost: assemblyAdjustedCost,
				synths:       nsynths,
				pcrs:         npcrs,
			}
			rlog.Debugf("Create filled assembly a[%d]; %v",
				selectedAssembliesStart+ai+1, filledAssembly)

			filled = append(filled, filledAssembly)
		}
	}
	return filled
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
