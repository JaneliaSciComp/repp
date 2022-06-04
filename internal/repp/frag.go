package repp

import (
	"fmt"
	"math"
	"os"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/jinzhu/copier"
)

var (
	// madePrimers, formerly made primers
	madePrimers = make(map[string][]Primer)

	// primerErrs, errors found during prior builds
	primerErrs = make(map[string]error)
)

// fragType is the Frag building type to be used in the assembly
type fragType int

const (
	// linear fragment, ie the type of a fragment as it was uploaded submitted and without PCR/synthesis
	linear fragType = iota

	// circular is a circular sequence of DNA, e.g.: many of Addgene's plasmids
	circular

	// pcr fragments are those prepared by pcr, often a subselection of their parent plasmid
	pcr

	// synthetic fragments are those that will be fully synthesized (eg: gBlocks)
	synthetic
)

// Frag is a single building block stretch of DNA for assembly
type Frag struct {
	// ID is a unique identifier for this fragment
	ID string `json:"id,omitempty"`

	// type of the fragment in string representation for export
	Type string `json:"type"`

	// Cost to make the fragment
	Cost float64 `json:"cost"`

	// fragment/plasmid's sequence
	Seq string `json:"seq,omitempty"`

	// sequence of a pcr fragment after PCR's addition of bp
	PCRSeq string `json:"pcrSeq,omitempty"`

	// primers necessary to create this (if pcr fragment)
	Primers []Primer `json:"primers,omitempty"`

	// fragType of this fragment. circular | pcr | synthetic | existing
	fragType fragType

	// uniqueID of a match, ID + the start index % seq-length
	// used identified that catches nodes that cross the zero-index
	uniqueID string

	// fullSeq is the entire seq of the Frag/fragment as it was read in (for forward engineering)
	fullSeq string

	// db that the frag came from
	db DB

	// start of this Frag on the target plasmid
	start int

	// end of this Frag on the target plasmid
	end int

	// start of the frag's first feature
	featureStart int

	// end of the frag's last covered feature
	featureEnd int

	// assemblies that span from this Frag to the end of the plasmid
	assemblies []assembly

	// build configuration
	conf *config.Config
}

// ranged: a stretch that a primer spans relative to the target sequence
type ranged struct {
	// start of the primer's range
	start int

	// end of the primer's range
	end int
}

// Primer is a single Primer used to create a PCR fragment
type Primer struct {
	// Seq of the primer (in 5' to 3' direction)
	Seq string `json:"seq"`

	// Strand of the primer; true if top strand, false if complement
	Strand bool `json:"strand"`

	// Penalty score
	Penalty float64 `json:"penalty"`

	// PairPenalty score from primer3
	PairPenalty float64 `json:"pairPenalty"`

	// Tm of the primer
	Tm float64 `json:"tm"`

	// GC % max
	GC float64 `json:"gc"`

	// Range that the primer spans on the fragment
	Range ranged `json:"-"`
}

// newFrag creates a Frag from a match
func newFrag(m match, conf *config.Config) *Frag {
	fType := pcr
	if m.circular {
		fType = circular
	}

	return &Frag{
		ID:       m.entry,
		uniqueID: m.uniqueID,
		Seq:      strings.ToUpper(m.seq),
		start:    m.queryStart,
		end:      m.queryEnd,
		db:       m.db,
		conf:     conf,
		fragType: fType,
	}
}

// newFlags is the plural of newFlag
func newFrags(matches []match, conf *config.Config) []*Frag {
	min := conf.FragmentsMinHomology
	max := conf.FragmentsMaxHomology

	var frags []*Frag
	for _, m := range matches {
		f := newFrag(m, conf)

		// try and shrink to avoid a duplicate junction with self
		selfJunction := f.junction(f, min, max)
		if selfJunction != "" {
			f.end -= len(selfJunction)
			if f.end-f.start < conf.PcrMinLength {
				continue
			}
			f.Seq = f.Seq[:len(f.Seq)-len(selfJunction)]
		}

		frags = append(frags, f)
	}
	return frags
}

// copy returns a deep dopy of a Frag. used because nodes are mutated
// during assembly filling, and we don't want primers being shared between
// nodes in different assemblies
func (f *Frag) copy() (newFrag *Frag) {
	newFrag = &Frag{}
	if err := copier.Copy(newFrag, f); err != nil {
		stderr.Fatal(err)
	}

	return
}

// cost returns the estimated cost of a fragment. Combination of source and preparation
func (f *Frag) cost(procure bool) (c float64) {
	if procure {
		c += f.db.Cost
	}

	if f.fragType == pcr && f.Primers != nil {
		// cost of primers plus the cost of a single PCR reaction
		c += float64(len(f.Primers[0].Seq)+len(f.Primers[1].Seq)) * f.conf.PcrBpCost
		c += f.conf.PcrRxnCost
	} else if f.fragType == synthetic {
		c += f.conf.SynthFragmentCost(len(f.Seq))
	}

	return
}

// distTo returns the distance between the start of this Frag and the end of the other.
// assumes that this Frag starts before the other
// will return a negative number if this Frag overlaps with the other and positive otherwise
func (f *Frag) distTo(other *Frag) (bpDist int) {
	return other.start - f.end
}

// overlapsViaPCR returns whether this Frag could overlap the other Frag through homology
// created via PCR
func (f *Frag) overlapsViaPCR(other *Frag) bool {
	return f.distTo(other) <= f.conf.PcrPrimerMaxEmbedLength
}

// overlapsViaHomology returns whether this Frag already has sufficient overlap with the
// other Frag without any preparation like PCR
func (f *Frag) overlapsViaHomology(other *Frag) bool {
	return f.distTo(other) <= -f.conf.FragmentsMinHomology
}

// synthDist returns the number of synthesized fragments that would need to be created
// between one Frag and another if the two were to be joined, with no existing
// fragments/nodes in-between, in an assembly
func (f *Frag) synthDist(other *Frag) (synthCount int) {
	dist := f.distTo(other)

	if f.overlapsViaPCR(other) {
		// if the dist is <MaxEmbedLength, we can PCR our way there
		// and add the mutated bp between the nodes with PCR
		return 0
	}

	floatDist := math.Max(1.0, float64(dist))

	// split up the distance between them by the max synthesized fragment size
	return int(math.Ceil(floatDist / float64(f.conf.SyntheticMaxLength)))
}

// costTo estimates the $ amount needed to get from this fragment
// to the other Frag passed, either by PCR or synthesis
//
// PCR is, of course, preferred if we can add homology between this and
// the other fragment with PCR homology arms alone
//
// Otherwise we find the total synthesis distance between this and
// the other fragment and divide that by the cost per bp of synthesized DNA
//
// This does not add in the cost of procurement which is in assembly.add()
func (f *Frag) costTo(other *Frag) (cost float64) {
	needsPCR := f.fragType == pcr || f.fragType == circular
	pcrNoHomology := 50.0 * f.conf.PcrBpCost // pcr no homology
	pcrHomology := (50.0 + float64(f.conf.FragmentsMinHomology)) * f.conf.PcrBpCost

	if other == f {
		if needsPCR {
			return pcrNoHomology
		}
		return 0
	}

	if f.overlapsViaPCR(other) {
		if f.overlapsViaHomology(other) {
			// there's already enough overlap between this Frag and the one being tested
			// estimating two primers, primer length assumed to be 25bp
			return pcrNoHomology
		}

		// we have to create some additional primer sequence to reach the next fragment
		// estimating here that we'll add half of minHomology to both sides
		return pcrHomology
	}

	// we need to create a new synthetic fragment to get from this fragment to the next
	// to account for both the bps between them as well as the additional bps we need to add
	// for homology between the two
	dist := f.distTo(other)
	dist += f.conf.FragmentsMinHomology * 2
	synthCost := f.conf.SynthFragmentCost(dist)

	// also account for whether this frag will require PCR
	if needsPCR {
		return synthCost + pcrNoHomology
	}
	return synthCost
}

// reach returns a slice of Frag indexes that overlap with, or are the first synth_count nodes
// away from this one within a slice of ordered nodes
func (f *Frag) reach(nodes []*Frag, i int, features bool) (reachable []int) {
	reachable = []int{}

	// accumulate the nodes that overlap with this one
	for {
		i++

		// we've run out of nodes
		if i >= len(nodes) {
			return reachable
		}

		if features && nodes[i].featureEnd <= f.featureEnd {
			continue
		} else if nodes[i].end < f.end {
			continue // fully engulfed
		}

		reachable = append(reachable, i)
		if nodes[i].uniqueID == f.uniqueID {
			break // do not go any further
		}
	}

	return
}

// junction checks for and returns any 100% identical homology between the end of this
// Frag and the start of the other. returns an empty string if there's no junction between them
func (f *Frag) junction(other *Frag, minHomology, maxHomology int) (junction string) {
	s1 := f.Seq
	if f.PCRSeq != "" {
		s1 = f.PCRSeq
	}

	s2 := other.Seq
	if other.PCRSeq != "" {
		s2 = other.PCRSeq
	}

	s1 = strings.ToUpper(s1)
	s2 = strings.ToUpper(s2)

	//      v-maxHomology from end    v-minHomology from end
	// ------------------------------------
	//                    -----------------------------
	start := len(s1) - maxHomology
	end := len(s1) - minHomology

	if start < 0 {
		start = 0
	}
	if end < 0 {
		end = 0
	}

	// for every possible start index
	for i := start; i <= end; i++ {
		// traverse from that index to the end of the seq
		for k, j := i, 0; k < len(s1); j, k = j+1, k+1 {
			if j >= len(s2) {
				break
			}

			if s1[k] != s2[j] {
				break
			}

			// we made it to the end of the sequence, there's a junction
			if k == len(s1)-1 {
				return s1[i:]
			}
		}
	}

	return
}

// synthTo returns synthetic fragments to get this Frag to the next.
// It creates a slice of building fragments that have homology against
// one another and are within the upper and lower synthesis bounds.
// target is the plasmid's full sequence. We need it to build up the target
// plasmid's sequence
func (f *Frag) synthTo(next *Frag, target string) (synths []*Frag) {
	jL := f.conf.FragmentsMinHomology // junction length

	// check whether we need to make synthetic fragments to get
	// to the next fragment in the assembly
	synCount := f.synthDist(next) // fragment count
	if synCount == 0 {
		return nil
	}

	tL := len(target)               // length of the full target plasmid
	fL := f.distTo(next) / synCount // each fragment's length
	fL += jL * 2                    // account for homology on either end of each synthetic fragment
	if f.conf.SyntheticMinLength > fL {
		// need to synthesize at least Synthesis.MinLength bps
		fL = f.conf.SyntheticMinLength
	}

	// add to self to account for sequence across the zero-index (when sequence subselecting)
	target = strings.ToUpper(target + target + target + target) // TODO remove this

	// slide along the range of sequence to create synthetic fragments
	// and create one at each point, each w/ jL for the fragment
	// before and after it
	synths = []*Frag{}
	start := f.end - jL + tL // start w/ homology, move left
	for len(synths) < synCount {
		end := start + fL + 1
		seq := target[start:end]

		// check for a hairpin in the junction and shift this fragment's synthesis
		// to the right if a hairpin is found
		for hairpin(seq[len(seq)-jL:], f.conf) > f.conf.FragmentsMaxHairpinMelt {
			end += jL / 2
			seq = target[start:end]
		}

		synths = append(synths, &Frag{
			ID:       fmt.Sprintf("%s-%s-synthesis-%d", f.ID, next.ID, len(synths)+1),
			Seq:      seq,
			start:    start,
			end:      end,
			fragType: synthetic,
			conf:     f.conf,
		})

		start = end - jL
	}

	return
}

// setPrimers creates primers against a Frag and returns an error if:
//	1. the primers have an unacceptably high primer3 penalty score
//	2. the primers have off-targets in their source plasmid/fragment
func (f *Frag) setPrimers(last, next *Frag, seq string, conf *config.Config) (err error) {
	pHash := primerHash(last, f, next)
	if oldPrimers, contained := madePrimers[pHash]; contained {
		f.Primers = oldPrimers
		mutatePrimers(f, seq, 0, 0) // set PCRSeq
		return nil
	}

	if oldErr, contained := primerErrs[pHash]; contained {
		return oldErr
	}

	psExec := newPrimer3(last, f, next, seq, conf)

	// make input file and write to the fs
	// find how many bp of additional sequence need to be added
	// to the left and right primers (too large for primer3_core)
	addLeft, addRight, err := psExec.input(
		conf.FragmentsMinHomology,
		conf.FragmentsMaxHomology,
		conf.PcrPrimerMaxEmbedLength,
		conf.PcrMinLength,
		conf.PcrBufferLength,
	)
	if err != nil {
		primerErrs[pHash] = err
		return
	}

	if err = psExec.run(); err != nil {
		primerErrs[pHash] = err
		return
	}

	if err = psExec.parse(seq); err != nil {
		primerErrs[pHash] = err
		return
	}

	// update Frag's range, and add additional bp to the left and right primer if it wasn't included in the primer3 output
	mutatePrimers(f, seq, addLeft, addRight)

	// make sure the fragment's length is still long enough for PCR
	if len(f.PCRSeq) < conf.PcrMinLength {
		err = fmt.Errorf(
			"failed to execute primer3: %s is %dbp, needs to be > %dbp",
			f.ID,
			f.end-f.start,
			conf.PcrMinLength,
		)
		f.Primers = nil
		primerErrs[pHash] = err
		return
	}

	// 1. check for whether the primers have too have a pair penalty score
	if f.Primers[0].PairPenalty > conf.PcrPrimerMaxPairPenalty {
		err = fmt.Errorf(
			"primers have pair primer3 penalty score of %f, should be less than %f:\f%+v\f%+v",
			f.Primers[0].PairPenalty,
			conf.PcrPrimerMaxPairPenalty,
			f.Primers[0],
			f.Primers[1],
		)
		f.Primers = nil
		primerErrs[pHash] = err
		return
	}

	// 2. check for whether either of the primers have an off-target/mismatch
	var mismatchExists bool
	var mm match

	if f.fullSeq != "" {
		// we have the full sequence (it was included in the forward design)
		mismatchResult := seqMismatch(f.Primers, f.ID, f.fullSeq, conf)
		mismatchExists = mismatchResult.wasMismatch
		mm = mismatchResult.m
		err = mismatchResult.err
	} else if f.db.Path != "" {
		// otherwise, query the fragment from the DB (try to find it) and then check for mismatches
		mismatchResult := parentMismatch(f.Primers, f.ID, f.db, conf)
		mismatchExists = mismatchResult.wasMismatch
		mm = mismatchResult.m
		err = mismatchResult.err
	}

	if err != nil {
		f.Primers = nil
		primerErrs[pHash] = err
		return err
	}
	if mismatchExists {
		err = fmt.Errorf(
			"found a mismatching sequence %s for primers: %s, %s",
			mm.seq,
			f.Primers[0].Seq,
			f.Primers[1].Seq,
		)
		f.Primers = nil
		primerErrs[pHash] = err
		return
	}

	f.fragType = pcr

	os.Remove(psExec.in.Name()) // delete the temporary input and output files
	os.Remove(psExec.out.Name())

	madePrimers[pHash] = f.Primers

	return
}

// mutatePrimers adds additional bp to the sides of a Frag
// if there was additional homology bearing sequence that we were unable
// to add through primer3 alone
//
// it also updates the range, start + end, of the Frag to match that of the primers
//
// returning Frag for testing
func mutatePrimers(f *Frag, seq string, addLeft, addRight int) *Frag {
	sl := len(seq)
	seq = strings.ToUpper(seq + seq + seq + seq) // TODO

	// change the Frag's start and end index to match those of the start and end index
	// of the primers, since the range may have shifted to get better primers
	f.start = f.Primers[0].Range.start
	f.end = f.Primers[1].Range.end

	// update fragment sequence
	f.Seq = seq[f.start+sl : f.end+sl+1]

	// add bp to the left/FWD primer to match the fragment to the left
	if addLeft > 0 {
		oldStart := f.Primers[0].Range.start + sl
		f.Primers[0].Seq = seq[oldStart-addLeft:oldStart] + f.Primers[0].Seq
		f.Primers[0].Range.start -= addLeft
	}

	// add bp to the right/REV primer to match the fragment to the right
	if addRight > 0 {
		oldEnd := f.Primers[1].Range.end + sl
		f.Primers[1].Seq = reverseComplement(seq[oldEnd+1:oldEnd+addRight+1]) + f.Primers[1].Seq
		f.Primers[1].Range.end += addRight
	}

	// update fragment sequence
	f.PCRSeq = seq[f.Primers[0].Range.start+sl : f.Primers[1].Range.end+sl+1]

	return f
}

// String returns a string representation of a fragment's type
func (t fragType) String() string {
	return []string{"linear", "plasmid", "pcr", "synthetic"}[t]
}

// fragsCost returns the total cost of a slice of frags. Just the summation of their costs
func fragsCost(frags []*Frag) (cost float64) {
	for _, f := range frags {
		cost += f.cost(true)
	}

	return
}

// primerHash returns a unique hash for a PCR run
func primerHash(last, f, next *Frag) string {
	return fmt.Sprintf("%s%d%d%d%d", f.uniqueID, last.end, f.start, f.end, next.start)
}
