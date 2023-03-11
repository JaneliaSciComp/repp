package repp

import (
	"fmt"
	"os"
	"regexp"
	"sort"
	"strings"
	"text/tabwriter"

	"github.com/Lattice-Automation/repp/internal/config"
)

// enzyme is a single enzyme that can be used to linearize a backbone before
// inserting a sequence.
type enzyme struct {
	name         string
	recog        string
	seqCutIndex  int
	compCutIndex int
}

// cut is a binding index and the length of the overhang after digestion
type cut struct {
	index  int
	strand bool
	enzyme enzyme
}

// Backbone is for information on a linearized backbone in the output payload
type Backbone struct {
	// URL of the backbone fragment's source
	URL string `json:"url"`

	// Seq is the sequence of the backbone (unlinearized)
	Seq string `json:"seq"`

	// Enzymes is the list of enzymes names used to linearize the backbone
	Enzymes []string `json:"enzymes"`

	// cutsites are the indexes where the backbone was cleaved
	Cutsites []int `json:"recognitionIndex"`

	// Strands of each cut direction. True if fwd, False if rev direction
	Strands []bool `json:"strands"`
}

// parses a recognition sequence into a hangInd, cutInd for overhang calculation.
func newEnzyme(name, recogSeq string) enzyme {
	cutIndex := strings.Index(recogSeq, "^")
	hangIndex := strings.Index(recogSeq, "_")

	if cutIndex < hangIndex {
		hangIndex--
	} else {
		cutIndex--
	}

	recogSeq = strings.Replace(recogSeq, "^", "", -1)
	recogSeq = strings.Replace(recogSeq, "_", "", -1)

	return enzyme{
		name:         name,
		recog:        recogSeq,
		seqCutIndex:  cutIndex,
		compCutIndex: hangIndex,
	}
}

// digest a Frag (backbone) with an enzyme's first recogition site
//
// remove the 5' end of the fragment post-cleaving. it will be degraded.
// keep exposed 3' ends. good visual explanation:
// https://warwick.ac.uk/study/csde/gsp/eportfolio/directory/pg/lsujcw/gibsonguide/
func digest(frag *Frag, enzymes []enzyme) (digested *Frag, backbone *Backbone, err error) {
	wrappedBp := 38 // largest current recognition site in the list of enzymes
	if len(frag.Seq) < wrappedBp {
		return &Frag{}, &Backbone{}, fmt.Errorf("%s is too short for digestion", frag.ID)
	}

	frag.Seq = strings.ToUpper(frag.Seq)
	firstHalf := frag.Seq[:len(frag.Seq)/2]
	secondHalf := frag.Seq[len(frag.Seq)/2:]
	if firstHalf == secondHalf {
		// it's a circular fragment that's doubled in the database
		frag.Seq = frag.Seq[:len(frag.Seq)/2] // undo the doubling of sequence for circular parts
	}

	// find all the cutsites
	cuts, lengths := cutsites(frag.Seq, enzymes)

	// none found
	if len(cuts) == 0 {
		enzymeNames := []string{}
		for _, enzyme := range enzymes {
			enzymeNames = append(enzymeNames, enzyme.name)
		}
		return &Frag{}, &Backbone{}, fmt.Errorf("no %s cutsites found in %s", strings.Join(enzymeNames, ","), frag.ID)
	}

	// only one cutsite
	if len(cuts) == 1 {
		cut := cuts[0]

		overhangLength := cut.enzyme.seqCutIndex - cut.enzyme.compCutIndex
		digestedSeq := ""

		if overhangLength >= 0 {
			cutIndex := (cut.index + cut.enzyme.seqCutIndex) % len(frag.Seq)
			digestedSeq = frag.Seq[cutIndex:] + frag.Seq[:cutIndex]
		} else {
			bottomIndex := (cut.index + cut.enzyme.seqCutIndex) % len(frag.Seq)
			topIndex := (cut.index + cut.enzyme.compCutIndex) % len(frag.Seq)
			digestedSeq = frag.Seq[topIndex:] + frag.Seq[:bottomIndex]
		}

		return &Frag{
				ID:       frag.ID,
				uniqueID: "backbone",
				Seq:      digestedSeq,
				fragType: linear,
				db:       frag.db,
			},
			&Backbone{
				Seq:      frag.Seq,
				Enzymes:  []string{cut.enzyme.name},
				Cutsites: []int{cut.index},
				Strands:  []bool{cut.strand},
			},
			nil
	}

	// find the largest band
	largestBand := 0
	for i, bandLength := range lengths {
		if bandLength > lengths[largestBand] {
			largestBand = i
		}
	}

	// find the enzyme from the start and end of the largest band
	cut1 := cuts[largestBand]
	cut2 := cuts[(largestBand+1)%len(lengths)]
	doubled := frag.Seq + frag.Seq

	cut1Index := cut1.index + cut1.enzyme.compCutIndex
	if cut1.enzyme.seqCutIndex-cut1.enzyme.compCutIndex < 0 {
		cut1Index = cut1.index + cut1.enzyme.seqCutIndex
	}

	cut2Index := cut2.index + cut2.enzyme.compCutIndex
	if cut2.enzyme.seqCutIndex-cut2.enzyme.compCutIndex < 0 {
		cut2Index = cut2.index + cut2.enzyme.seqCutIndex
	}

	if cut2Index < cut1Index {
		cut2Index += len(frag.Seq)
	}

	digestedSeq := doubled[cut1Index:cut2Index]

	return &Frag{
			ID:       frag.ID,
			uniqueID: "backbone",
			Seq:      digestedSeq,
			fragType: linear,
			db:       frag.db,
		},
		&Backbone{
			Seq:      frag.Seq,
			Enzymes:  []string{cut1.enzyme.name, cut2.enzyme.name},
			Cutsites: []int{cut1Index, cut2Index},
			Strands:  []bool{cut1.strand, cut2.strand},
		},
		nil
}

// cutsites finds all the cutsites of a list of enzymes against a target sequence
// also returns the lengths of each "band" of DNA after digestion. Each band length
// corresponds to the band formed with the start of the enzyme at the same index in cuts
func cutsites(seq string, enzymes []enzyme) (cuts []cut, lengths []int) {
	s := seq
	rcs := reverseComplement(seq)

	for _, enzyme := range enzymes {
		regexRecognition := recogRegex(enzyme.recog)
		reg := regexp.MustCompile(regexRecognition)

		for _, submatch := range reg.FindAllStringSubmatchIndex(s, -1) {
			index := submatch[0]
			cuts = append(cuts, cut{index: index, enzyme: enzyme, strand: true})
		}

		// if it's a palindrome enzyme, don't scan over it again
		if reverseComplement(regexRecognition) == regexRecognition {
			continue
		}

		for _, submatch := range reg.FindAllStringSubmatchIndex(rcs, -1) {
			index := submatch[0]
			index = len(seq) - index - len(enzyme.recog)
			cuts = append(cuts, cut{index: index, enzyme: enzyme, strand: false})
		}
	}

	sort.Slice(cuts, func(i, j int) bool {
		return cuts[i].index < cuts[j].index
	})

	for i, c := range cuts {
		next := (i + 1) % len(cuts)
		bandLength := (cuts[next].index - c.index + len(seq)) % len(seq)
		lengths = append(lengths, bandLength)
	}

	return
}

// recogRegex turns a recognition sequence into a regex sequence for searching
// sequence for searching the sequence for digestion sites.
func recogRegex(recog string) (decoded string) {
	regexDecode := map[rune]string{
		'A': "A",
		'C': "C",
		'G': "G",
		'T': "T",
		'M': "(A|C)",
		'R': "(A|G)",
		'W': "(A|T)",
		'Y': "(C|T)",
		'S': "(C|G)",
		'K': "(G|T)",
		'H': "(A|C|T)",
		'D': "(A|G|T)",
		'V': "(A|C|G)",
		'B': "(C|G|T)",
		'N': "(A|C|G|T)",
		'X': "(A|C|G|T)",
	}

	var regexDecoder strings.Builder
	for _, c := range recog {
		regexDecoder.WriteString(regexDecode[c])
	}

	return regexDecoder.String()
}

// NewEnzymeDB returns a new copy of the enzymes db.
func NewEnzymeDB() *kv {
	return newKV(config.EnzymeDB)
}

// PrintEnzymes writes enzymes that are similar in queried name to stdout.
// if multiple enzyme names include the enzyme name, they are all returned.
// otherwise a list of enzyme names are returned (those beneath a levenshtein distance cutoff).
func PrintEnzymes(enzyme string) {
	f := NewEnzymeDB()

	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)

	if enzyme == "" {
		enzymeNames := make([]string, len(f.contents))
		i := 0
		for name := range f.contents {
			enzymeNames[i] = name
			i++
		}
		sort.Strings(enzymeNames)

		for _, name := range enzymeNames {
			fmt.Fprintf(w, "%s\t%s\n", name, f.contents[name])
		}
		w.Flush()
		return
	}

	// if there's an exact match, just log that one
	if seq, exists := f.contents[enzyme]; exists {
		fmt.Printf("%s	%s\n", enzyme, seq)
		return
	}

	ldCutoff := 2
	containing := []string{}
	lowDistance := []string{}

	for fName, fSeq := range f.contents {
		if strings.Contains(fName, enzyme) {
			containing = append(containing, fName+"\t"+fSeq)
		} else if len(fName) > ldCutoff && ld(enzyme, fName, true) <= ldCutoff {
			lowDistance = append(lowDistance, fName+"\t"+fSeq)
		}
	}

	if len(containing) < 3 {
		lowDistance = append(lowDistance, containing...)
		containing = []string{} // clear
	}
	if len(containing) > 0 {
		fmt.Fprint(w, strings.Join(containing, "\n"))
	} else if len(lowDistance) > 0 {
		fmt.Fprint(w, strings.Join(lowDistance, "\n"))
	} else {
		fmt.Fprintf(w, "failed to find any enzymes for %s", enzyme)
	}
	if _, err := w.Write([]byte("\n")); err != nil {
		rlog.Fatal(err)
	}
	w.Flush()
}

// AddEnzymes the enzyme's seq in the database (or create if it isn't in the enzyme db).
func AddEnzymes(name, inputSeq string) {
	f := NewEnzymeDB()

	invalidChars := regexp.MustCompile(`[^ATGCMRWYSKHDVBNX_\^]`)
	seq := invalidChars.ReplaceAllString(strings.ToUpper(inputSeq), "")

	if strings.Count(seq, "^") != 1 || strings.Count(seq, "_") != 1 {
		rlog.Fatal("%s is not a valid enzyme recognition sequence. see 'repp find enzyme --help'\n", seq)
	}

	f.contents[name] = seq
	if err := f.save(); err != nil {
		rlog.Fatal(err)
	}
}

// DeleteEnzyme deletes one or more enzymes from the database
func DeleteEnzyme(enzyme string) (bresult bool, err error) {
	f := NewEnzymeDB()

	if _, contained := f.contents[enzyme]; !contained {
		return false, fmt.Errorf("failed to find %s in the enzymes database\n", enzyme)
	}

	delete(f.contents, enzyme)
	if err := f.save(); err != nil {
		return false, err
	}
	bresult = true

	return
}

func getValidEnzymes(enzymeNames []string) (enzymes []enzyme, err error) {
	enzymeDB := NewEnzymeDB()
	for _, enzymeName := range enzymeNames {
		if cutseq, exists := enzymeDB.contents[enzymeName]; exists {
			enzymes = append(enzymes, newEnzyme(enzymeName, cutseq))
		} else {
			return enzymes, fmt.Errorf(
				`failed to find enzyme with name %s use "repp enzymes" for a list of recognized enzymes`,
				enzymeName,
			)
		}
	}

	return
}
