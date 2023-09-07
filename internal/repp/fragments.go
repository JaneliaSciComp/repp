package repp

import (
	"fmt"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
)

// PrintFragment logs the building fragment with the name passed.
func PrintFragment(name string, dbNames []string) {
	dbs, err := getRegisteredDBs(dbNames)
	if err != nil {
		rlog.Fatal(err)
	}
	frag, err := queryDatabases(name, dbs)
	if err != nil {
		rlog.Fatal(err)
	}
	if frag.fragType == circular {
		frag.Seq = frag.Seq[:len(frag.Seq)/2]
	}
	fmt.Printf("%s\t%s\n%s\n", name, frag.db.Name, frag.Seq)
}

// AssembleFragments assembles a list of building fragments in order
func AssembleFragments(assemblyParams AssemblyParams, conf *config.Config) {

	// read in the constituent fragments
	frags, err := read(assemblyParams.GetIn(), false, false)
	if err != nil {
		rlog.Fatal(err)
	}
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
	// add in the backbone if it was provided
	if backboneFrag.ID != "" {
		frags = append([]*Frag{backboneFrag}, frags...)
	}

	// set the conf property on each frag
	for _, f := range frags {
		f.conf = conf
	}

	target, solution := fragments(frags, conf)

	primersDB := readOligos(assemblyParams.GetPrimersDBLocations(), primerIDPrefix, false)
	synthFragsDB := readOligos(assemblyParams.GetSynthFragsDBLocations(), synthFragIDPrefix, true)

	// write the single list of fragments as a possible solution to the output file
	if _, err := writeResult(
		assemblyParams.GetOut(),
		assemblyParams.GetOutputFormat(),
		assemblyParams.GetIn(),
		target.Seq,
		[][]*Frag{solution},
		primersDB,
		synthFragsDB,
		len(target.Seq),
		0,
		backboneMeta,
		conf,
	); err != nil {
		rlog.Fatal(err)
	}
}

// fragments pieces together a list of fragments into a single plasmid
// with the fragments in the order and orientation specified
func fragments(frags []*Frag, conf *config.Config) (target *Frag, solution []*Frag) {
	// piece together the adjacent fragments
	if len(frags) < 1 {
		rlog.Fatal("failed: no fragments to assemble")
	}

	// anneal the fragments together, shift their junctions and create the plasmid sequence
	vecSeq := annealFragments(conf.FragmentsMinHomology, conf.FragmentsMaxHomology, frags)

	// create the assumed target plasmid object
	target = &Frag{
		Seq:      vecSeq,
		fragType: circular,
	}

	// create an assembly out of the frags (to fill/convert to fragments with primers)
	a := assembly{frags: frags}
	solution, err := a.fill(target.Seq, conf)
	if err != nil {
		rlog.Fatal(err)
	}

	return target, solution
}

// annealFragments shifts the start and end of junctions that overlap one another
func annealFragments(min, max int, frags []*Frag) (vec string) {
	// set the start, end, and plasmid sequence
	// add all of each frags seq to the plasmid sequence, minus the region overlapping the next
	var vecSeq strings.Builder
	for i, f := range frags {
		next := frags[(i+1)%len(frags)]
		// if we're on the last fragment, mock the first one further along the plasmid
		if i == len(frags)-1 {
			nextSeq := next.Seq
			if next.PCRSeq != "" {
				nextSeq = next.PCRSeq
			}
			next = &Frag{
				Seq:   nextSeq,
				start: next.start + vecSeq.Len(),
				end:   next.end + vecSeq.Len(),
			}
		}

		j := len(f.junction(next, min, max)) // junction length

		fragSeq := f.Seq
		if f.PCRSeq != "" {
			fragSeq = f.PCRSeq
		}
		contrib := fragSeq[0 : len(fragSeq)-j] // frag's contribution to plasmid

		// correct for this Frag's overlap with the next Frag
		f.start = vecSeq.Len()
		f.end = f.start + len(fragSeq) - 1

		// add this Frag's sequence onto the accumulated plasmid sequence
		vecSeq.WriteString(contrib)
	}

	return vecSeq.String()
}

// validateJunctions checks each fragment and confirms that it has sufficient homology
// with its adjacent fragments and that the match is exact. Largely for testing
func validateJunctions(frags []*Frag, conf *config.Config) error {
	for i, f := range frags {
		next := frags[(i+1)%len(frags)]
		j := f.junction(next, conf.FragmentsMinHomology, conf.FragmentsMaxHomology+1)
		if j == "" {
			s1 := f.getFragSeq()
			s2 := next.getFragSeq()

			currID := f.ID
			nextID := next.ID
			return fmt.Errorf("no junction found between %s and %s\n%s\n\n%s", currID, nextID, s1, s2)
		}
	}

	return nil
}
