package repp

import (
	"fmt"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/spf13/cobra"
)

// FragmentListCmd logs the building fragment with the name passed.
func FragmentListCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		cmd.Help()
		stderr.Fatalln("\nno fragment name passed.")
	}
	name := args[0]

	flags, _ := parseCmdFlags(cmd, args, false)
	frag, err := queryDatabases(name, flags.dbs)
	if err != nil {
		stderr.Fatalln(err)
	}
	if frag.fragType == circular {
		frag.Seq = frag.Seq[:len(frag.Seq)/2]
	}

	fmt.Printf("%s\t%s\n%s\n", name, frag.db.GetName(), frag.Seq)
}

// FragmentsCmd accepts a cobra commands and assembles a list of building fragments in order
func FragmentsCmd(cmd *cobra.Command, args []string) {
	flags, conf := parseCmdFlags(cmd, args, true)

	// read in the constituent fragments
	frags, err := read(flags.in, false)
	if err != nil {
		stderr.Fatalln(err)
	}

	// add in the backbone if it was provided
	if flags.backbone.ID != "" {
		frags = append([]*Frag{flags.backbone}, frags...)
	}

	// set the conf property on each frag
	for _, f := range frags {
		f.conf = conf
	}

	target, solution := fragments(frags, conf)

	// write the single list of fragments as a possible solution to the output file
	writeJSON(
		flags.out,
		flags.in,
		target.Seq,
		[][]*Frag{solution},
		len(target.Seq),
		0,
		flags.backboneMeta,
		conf,
	)
}

// fragments pieces together a list of fragments into a single plasmid
// with the fragments in the order and orientation specified
func fragments(frags []*Frag, conf *config.Config) (target *Frag, solution []*Frag) {
	// piece together the adjacent fragments
	if len(frags) < 1 {
		stderr.Fatalln("failed: no fragments to assemble")
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
		stderr.Fatalln(err)
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
			s1 := f.Seq
			if f.PCRSeq != "" {
				s1 = f.PCRSeq
			}

			s2 := next.Seq
			if next.PCRSeq != "" {
				s2 = next.PCRSeq
			}

			left := f.ID
			right := next.ID
			return fmt.Errorf("no junction found between %s and %s\n%s\n\n%s", left, right, s1, s2)
		}
	}

	return nil
}
