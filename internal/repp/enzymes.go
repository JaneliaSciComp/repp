package repp

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"os"
	"regexp"
	"sort"
	"strings"
	"text/tabwriter"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/spf13/cobra"
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
				URL:      parseURL(frag.ID, frag.db),
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
			URL:      parseURL(frag.ID, frag.db),
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

// EnzymeDB is a struct for accessing repps enzymes db.
type EnzymeDB struct {
	// enzymes is a map between a enzymes name and its sequence
	enzymes map[string]string
}

// NewEnzymeDB returns a new copy of the enzymes db.
func NewEnzymeDB() *EnzymeDB {
	enzymeFile, err := os.Open(config.EnzymeDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	scanner := bufio.NewScanner(enzymeFile)
	enzymes := make(map[string]string)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		enzymes[columns[0]] = columns[1] // enzyme name = enzyme seq
	}

	if err := enzymeFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	return &EnzymeDB{enzymes: enzymes}
}

// ReadCmd returns enzymes that are similar in name to the enzyme name requested.
// if multiple enzyme names include the enzyme name, they are all returned.
// otherwise a list of enzyme names are returned (those beneath a levenshtein distance cutoff).
func (f *EnzymeDB) ReadCmd(cmd *cobra.Command, args []string) {
	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)

	if len(args) < 1 {
		enzymeNames := make([]string, len(f.enzymes), len(f.enzymes))
		i := 0
		for name := range f.enzymes {
			enzymeNames[i] = name
			i++
		}
		sort.Strings(enzymeNames)

		for _, name := range enzymeNames {
			fmt.Fprintf(w, "%s\t%s\n", name, f.enzymes[name])
		}
		w.Flush()
		return
	}

	name := args[0]

	// if there's an exact match, just log that one
	if seq, exists := f.enzymes[name]; exists {
		fmt.Printf("%s	%s\n", name, seq)
		return
	}

	ldCutoff := 2
	containing := []string{}
	lowDistance := []string{}

	for fName, fSeq := range f.enzymes {
		if strings.Contains(fName, name) {
			containing = append(containing, fName+"\t"+fSeq)
		} else if len(fName) > ldCutoff && ld(name, fName, true) <= ldCutoff {
			lowDistance = append(lowDistance, fName+"\t"+fSeq)
		}
	}

	if len(containing) < 3 {
		lowDistance = append(lowDistance, containing...)
		containing = []string{} // clear
	}
	if len(containing) > 0 {
		fmt.Fprintf(w, strings.Join(containing, "\n"))
	} else if len(lowDistance) > 0 {
		fmt.Fprintf(w, strings.Join(lowDistance, "\n"))
	} else {
		fmt.Fprintf(w, fmt.Sprintf("failed to find any enzymes for %s", name))
	}
	w.Write([]byte("\n"))
	w.Flush()
}

// SetCmd the enzyme's seq in the database (or create if it isn't in the enzyme db).
func (f *EnzymeDB) SetCmd(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		cmd.Help()
		stderr.Fatalln("expecting two args: a name and recognition sequence.")
	}

	name := args[0]
	seq := args[1]
	if len(args) > 2 {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]
	}
	seq = strings.ToUpper(seq)

	invalidChars := regexp.MustCompile("[^ATGCMRWYSKHDVBNX_\\^]")
	seq = invalidChars.ReplaceAllString(seq, "")

	if strings.Count(seq, "^") != 1 || strings.Count(seq, "_") != 1 {
		stderr.Fatalf("%s is not a valid enzyme recognition sequence. see 'repp find enzyme --help'\n", seq)
	}

	enzymeFile, err := os.Open(config.EnzymeDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	var output strings.Builder
	updated := false
	scanner := bufio.NewScanner(enzymeFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		if columns[0] == name {
			output.WriteString(fmt.Sprintf("%s	%s\n", name, seq))
			updated = true
		} else {
			output.WriteString(scanner.Text())
		}
	}

	// create from nothing
	if !updated {
		output.WriteString(fmt.Sprintf("%s	%s\n", name, seq))
	}

	if err := enzymeFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	if err := ioutil.WriteFile(config.EnzymeDB, []byte(output.String()), 0644); err != nil {
		stderr.Fatal(err)
	}

	if updated {
		fmt.Printf("updated %s in the enzymes database\n", name)
	}

	// update in memory
	f.enzymes[name] = seq
}

// DeleteCmd the enzyme from the database
func (f *EnzymeDB) DeleteCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		cmd.Help()
		stderr.Fatalf("\nexpecting an enzymes name.")
	}

	name := args[0]
	if len(args) > 1 {
		name = strings.Join(args, " ")
	}

	if _, contained := f.enzymes[name]; !contained {
		fmt.Printf("failed to find %s in the enzymes database\n", name)
	}

	enzymeFile, err := os.Open(config.EnzymeDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	var output strings.Builder
	deleted := false
	scanner := bufio.NewScanner(enzymeFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		if columns[0] != name {
			output.WriteString(scanner.Text())
		} else {
			deleted = true
		}
	}

	if err := enzymeFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	if err := ioutil.WriteFile(config.EnzymeDB, []byte(output.String()), 0644); err != nil {
		stderr.Fatal(err)
	}

	// delete from memory
	delete(f.enzymes, name)

	if deleted {
		fmt.Printf("deleted %s from the enzymes database\n", name)
	} else {
		fmt.Printf("failed to find %s in the enzymes database\n", name)
	}
}
