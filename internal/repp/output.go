package repp

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/Lattice-Automation/repp/internal/config"
	"go.uber.org/multierr"
)

// Solution is a single solution to build up the target plasmid.
type Solution struct {
	// Count is the number of fragments in this solution
	Count int `json:"count"`

	// Cost estimated from the primer and sequence lengths
	Cost float64 `json:"cost"`

	// Adjusted cost for synthentic fragments
	AdjustedCost float64 `json:"adjustedCost"`

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`

	// number of PCR fragments
	pcrFragsCount int

	// number of synthentic fragments
	synthFragsCount int
}

// Output is a struct containing design results for the assembly.
type Output struct {
	// Target's name. In >example_CDS FASTA its "example_CDS"
	Target string `json:"target"`

	// Target's sequence
	TargetSeq string `json:"seq"`

	// Time, ex: "2018-01-01 20:41:00"
	Time string `json:"time"`

	// Execution is the number of seconds it took to execute the command
	Execution float64 `json:"execution"`

	// Solutions builds
	Solutions []Solution `json:"solutions"`

	// Backbone is the user linearized a backbone fragment
	Backbone *Backbone `json:"backbone,omitempty"`
}

// writeResult
func writeResult(
	filename,
	format,
	targetName,
	targetSeq string,
	assemblies [][]*Frag,
	primersDB, synthFragsDB *oligosDB,
	backbone *Backbone,
	seconds float64,
	conf *config.Config,
) (*Output, error) {
	out, err := prepareSolutionsOutput(
		targetName,
		targetSeq,
		assemblies,
		backbone,
		seconds,
		conf,
	)
	if err != nil {
		return nil, err
	}
	if format == "CSV" {
		err = writeCSV(filename, fragmentBase(filename), primersDB, synthFragsDB, conf.IncludeFragLocationInStrategyOutput, out)
	} else {
		err = writeJSON(filename, out)
	}
	return out, err
}

// prepareSolutionsOutput turns a list of solutions into a Solution object.
func prepareSolutionsOutput(
	targetName,
	targetSeq string,
	assemblies [][]*Frag,
	backbone *Backbone,
	seconds float64,
	conf *config.Config,
) (out *Output, err error) {
	// store save time, using same format as log.Println https://golang.org/pkg/log/#Println
	t := time.Now() // https://gobyexample.com/time-formatting-parsing
	time := fmt.Sprintf(
		"%d/%02d/%02d %02d:%02d:%02d",
		t.Year(), t.Month(), t.Day(), t.Hour(), t.Minute(), t.Second(),
	)
	roundCost := func(cost float64) (float64, error) {
		return strconv.ParseFloat(fmt.Sprintf("%.2f", cost), 64)
	}

	// calculate final cost of the assembly and fragment count
	solutions := []Solution{}
	for _, assembly := range assemblies {
		assemblyCost := 0.0
		assemblyAdjustedCost := 0.0
		assemblyFragmentIDs := make(map[string]bool)
		gibson := false // whether it will be assembled via Gibson assembly
		hasPCR := false // whether there will be a batch PCR
		npcrs := 0
		nsynths := 0
		for _, f := range assembly {
			var fragCost, fragAdjustedCost float64
			if f.fragType != linear && f.fragType != circular {
				gibson = true
			}

			if f.fragType == pcr {
				hasPCR = true
				npcrs++
			} else {
				nsynths++
			}
			f.Type = f.fragType.String() // freeze fragment type

			// if it's already in the assembly, don't count cost twice
			if _, contained := assemblyFragmentIDs[f.ID]; f.ID != "" && contained {
				fragCost, fragAdjustedCost = f.cost(false)
			} else {
				fragCost, fragAdjustedCost = f.cost(true) // do not include procurement costs twice
				assemblyFragmentIDs[f.ID] = true
			}
			// round to two decimal places
			if f.Cost, err = roundCost(fragCost); err != nil {
				return nil, err
			}
			if f.AdjustedCost, err = roundCost(fragAdjustedCost); err != nil {
				return nil, err
			}

			// accumulate assembly cost
			assemblyCost += f.Cost
			assemblyAdjustedCost += f.AdjustedCost
		}

		if gibson {
			assemblyCost += conf.GibsonAssemblyCost + conf.GibsonAssemblyTimeCost
			assemblyAdjustedCost += conf.GibsonAssemblyCost + conf.GibsonAssemblyTimeCost
		}

		if hasPCR {
			assemblyCost += conf.PcrTimeCost
			assemblyAdjustedCost += conf.PcrTimeCost
		}

		solutionCost, err := roundCost(assemblyCost)
		if err != nil {
			return nil, err
		}
		solutionAdjustedCost, err := roundCost(assemblyAdjustedCost)
		if err != nil {
			return nil, err
		}

		solutions = append(solutions, Solution{
			Count:           len(assembly),
			Cost:            solutionCost,
			AdjustedCost:    solutionAdjustedCost,
			Fragments:       assembly,
			pcrFragsCount:   npcrs,
			synthFragsCount: nsynths,
		})
	}

	// sort solutions in increasing fragment count order
	sort.Slice(solutions, func(i, j int) bool {
		return solutions[i].Count < solutions[j].Count
	})

	if backbone.Seq == "" {
		backbone = nil
	}

	out = &Output{
		Time:      time,
		Target:    targetName,
		TargetSeq: strings.ToUpper(targetSeq),
		Execution: seconds,
		Solutions: solutions,
		Backbone:  backbone,
	}

	return out, nil
}

// writeCSV writes solutions as csv.
// The results are output to two csv files;
// one containing the strategy and the other one the reagents
func writeCSV(filename, fragmentIDBase string,
	existingPrimers, existingSynthFrags *oligosDB,
	withFragLocation bool,
	out *Output) (err error) {

	reagentsFilename := resultFilename(filename, "reagents")
	strategyFilename := resultFilename(filename, "strategy")

	reagentsFile, err := os.Create(reagentsFilename)
	if err != nil {
		return err
	}
	defer reagentsFile.Close()

	strategyFile, err := os.Create(strategyFilename)
	if err != nil {
		return err
	}
	defer strategyFile.Close()

	strategyCSVWriter := csv.NewWriter(strategyFile)
	// write timestamp
	_, err = fmt.Fprintf(strategyFile, "# %s\n", out.Time)
	if err != nil {
		return err
	}

	reagentsCSVWriter := csv.NewWriter(reagentsFile)
	// Write the strategy headers
	var headers []string
	if withFragLocation {
		headers = []string{
			"Frag ID",
			"Fwd Primer",
			"Rev Primer",
			"Template",
			"Size",
			"Match Pct",
			"Frag Start",
			"Frag End",
			"Template Start",
			"Template End",
			"GC%",
			"50 low GC%",
			"50 high GC%",
			"Homopolymer",
		}
	} else {
		headers = []string{
			"Frag ID",
			"Fwd Primer",
			"Rev Primer",
			"Template",
			"Size",
			"Match Pct",
			"GC%",
			"50 low GC%",
			"50 high GC%",
			"Homopolymer",
		}
	}
	err = strategyCSVWriter.Write(headers)
	if err != nil {
		return nil
	}
	// Write the reagents headers
	err = reagentsCSVWriter.Write([]string{
		"Reagent ID",
		"Seq",
		"Priming Region",
		"Tm",
		"Notes",
	})
	for si, s := range out.Solutions {
		snumber := si + 1
		// Write the solution cost and the number of fragments
		if _, err = fmt.Fprintf(strategyFile,
			"# Solution %d\n# Fragments:%d (%d - pcr, %d - synth)\n# Cost: %f, Adjusted Cost: %f\n",
			snumber,
			s.Count, s.pcrFragsCount, s.synthFragsCount,
			s.Cost, s.AdjustedCost); err != nil {
			return err
		}
		if _, err = fmt.Fprintf(reagentsFile, "# Solution %d\n", snumber); err != nil {
			return err
		}
		reagents := []oligo{}
		var newPrimerIndex int = 0
		var newSynthFragIndex int = 0

		newPrimers := newOligosDB(primerIDPrefix, false)
		newSynthFrags := newOligosDB(synthFragIDPrefix, true)

		var updatedPrimerDBs []*oligosDB = []*oligosDB{
			existingPrimers,
			newPrimers,
		}

		var updatedSynthFragsDBs []*oligosDB = []*oligosDB{
			existingSynthFrags,
			newSynthFrags,
		}

		for fi, f := range s.Fragments {
			fnumber := fi + 1
			var fwdPrimer, revPrimer Primer
			var synthSeq string

			fID := fmt.Sprintf("%s_%d_%s", fragmentIDBase, fnumber, fragTypeAsString(f.fragType))
			fwdPrimer, revPrimer = f.getPrimers()
			if fwdPrimer.Seq == "" && revPrimer.Seq == "" {
				synthSeq = f.Seq
			}

			fwdOligo := searchOligoDBs(fwdPrimer.Seq, updatedPrimerDBs)
			if !fwdOligo.isEmpty() {
				if !fwdOligo.hasID() {
					fwdOligo.assignNewOligoID(existingPrimers.getNewOligoID(newPrimerIndex))
					newPrimers.addOligo(fwdOligo)
					newPrimerIndex++
				}
				fwdOligo.primingRegion = fwdPrimer.PrimingRegion
				fwdOligo.tm = fwdPrimer.Tm
				fwdOligo.notes = fwdPrimer.Notes
				reagents = append(reagents, fwdOligo)
			}
			revOligo := searchOligoDBs(revPrimer.Seq, updatedPrimerDBs)
			if !revOligo.isEmpty() {
				if !revOligo.hasID() {
					revOligo.assignNewOligoID(existingPrimers.getNewOligoID(newPrimerIndex))
					newPrimers.addOligo(revOligo)
					newPrimerIndex++
				}
				revOligo.primingRegion = revPrimer.PrimingRegion
				revOligo.tm = revPrimer.Tm
				revOligo.notes = revPrimer.Notes
				reagents = append(reagents, revOligo)
			}
			var templateID string
			var matchRatio string
			var pcrSeqSize int
			var fragStart, fragEnd, templateStart, templateEnd string
			var gcContentCol string
			var min50GCContentCol string
			var max50GCContentCol string
			var homopolymerCol string
			if f.fragType == synthetic {
				synthReagent := searchOligoDBs(synthSeq, updatedSynthFragsDBs)
				if !synthReagent.hasID() {
					synthReagent.assignNewOligoID(existingSynthFrags.getNewOligoID(newSynthFragIndex))
					synthReagent.synth = true
					newSynthFrags.addOligo(synthReagent)
					newSynthFragIndex++
				}
				fID = synthReagent.id
				templateID = "N/A"
				matchRatio = "N/A"
				pcrSeqSize = len(f.Seq)
				fragStart = fmt.Sprintf("%d", f.start)
				fragEnd = fmt.Sprintf("%d", f.end)
				templateStart = "N/A"
				templateEnd = "N/A"
				reagents = append(reagents, synthReagent)
				synthFragScores := fragSeqQualityChecks(f.Seq)
				gcContentCol = fmt.Sprintf("%3.1f", synthFragScores.gcContent*100)
				min50GCContentCol = fmt.Sprintf("%3.1f", synthFragScores.min50WindowGCContent*100)
				max50GCContentCol = fmt.Sprintf("%3.1f", synthFragScores.max50WindowGCContent*100)
				homopolymerCol = strconv.Itoa(synthFragScores.longestHomopolymer)
			} else {
				templateID = fragmentBase(f.ID)
				matchRatio = fmt.Sprintf("%d", int(f.matchRatio*100))
				// for PCR fragments display the length including the overhanging primers
				pcrSeqSize = len(f.PCRSeq)
				if f.revCompFlag {
					fragStart = fmt.Sprintf("%d", f.end)
					if f.start >= len(out.TargetSeq) {
						fragEnd = fmt.Sprintf("%d(-)", f.start-len(out.TargetSeq))
					} else {
						fragEnd = fmt.Sprintf("%d", f.start)
					}
				} else {
					fragStart = fmt.Sprintf("%d", f.start)
					if f.end >= len(out.TargetSeq) {
						fragEnd = fmt.Sprintf("%d(+)", f.end-len(out.TargetSeq))
					} else {
						fragEnd = fmt.Sprintf("%d", f.end)
					}
				}
				if f.revCompTemplateFlag {
					templateStart = fmt.Sprintf("%d", f.templateEnd)
					templateEnd = fmt.Sprintf("%d", f.templateStart)
				} else {
					templateStart = fmt.Sprintf("%d", f.templateStart)
					templateEnd = fmt.Sprintf("%d", f.templateEnd)
				}
				gcContentCol = "N/A"
				min50GCContentCol = "N/A"
				max50GCContentCol = "N/A"
				homopolymerCol = "N/A"
			}
			fieldMapping := map[string]string{
				"Frag ID":        fID,
				"Fwd Primer":     fwdOligo.getIDOrDefault(false, "N/A"), // fwd primer
				"Rev Primer":     revOligo.getIDOrDefault(false, "N/A"), // rev primer
				"Template":       templateID,                            // template
				"Size":           strconv.Itoa(pcrSeqSize),
				"Match Pct":      matchRatio,
				"Frag Start":     fragStart,
				"Frag End":       fragEnd,
				"Template Start": templateStart,
				"Template End":   templateEnd,
				"GC%":            gcContentCol,
				"50 low GC%":     min50GCContentCol,
				"50 high GC%":    max50GCContentCol,
				"Homopolymer":    homopolymerCol,
			}
			var fields []string
			for _, h := range headers {
				fields = append(fields, fieldMapping[h])
			}
			if err = strategyCSVWriter.Write(fields); err != nil {
				return nil
			}
		}
		strategyCSVWriter.Flush()
		sort.Sort(sortedOligosByID(reagents))
		for _, r := range reagents {
			err = writeReagent(reagentsCSVWriter, r)
			if err != nil {
				rlog.Errorf("Error writing reagent %s: %v", r.id, err)
			}
		}
		reagentsCSVWriter.Flush()
	}

	return nil
}

func fragmentBase(filename string) string {
	baseNameFromFilename := fragIDComponents(filepath.Base(filename))[0]
	if len(baseNameFromFilename) > 10 {
		// truncate if name is too long
		return baseNameFromFilename[:10]
	} else {
		return baseNameFromFilename
	}
}

func fragIDComponents(fragID string) []string {
	var fragmentIDSeparators = " ,_-.;:"

	splitFunc := func(c rune) bool {
		return strings.ContainsRune(fragmentIDSeparators, c)
	}

	return strings.FieldsFunc(fragID, splitFunc)
}

func resultFilename(template, suffix string) string {
	ext := filepath.Ext(template)
	noExt := template[0 : len(template)-len(ext)]
	return noExt + "-" + suffix + ext
}

func writeReagent(csvWriter *csv.Writer, reagent oligo) (err error) {
	reagentID := reagent.getIDOrDefault(!reagent.isNew, "N/A") // mark the ID if this reagent already existed in the original manifest
	if reagentID != "" {
		var primingRegion, tm string
		if reagent.primingRegion == "" {
			primingRegion = "N/A"
			tm = "N/A"
		} else {
			primingRegion = reagent.primingRegion
			tm = fmt.Sprintf("%.2f", reagent.tm)
		}
		err = csvWriter.Write([]string{
			reagentID,
			reagent.seq,
			primingRegion,
			tm,
			reagent.notes,
		})
	}
	return
}

// writeJSON writes solutions as json.
func writeJSON(filename string, out *Output) (err error) {
	output, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		return fmt.Errorf("failed to serialize output: %v", err)
	}

	if err = os.WriteFile(filename, output, 0666); err != nil {
		return fmt.Errorf("failed to write the output: %v", err)
	}

	return
}

// writeFragsToFastaFile writes a slice of fragments to a FASTA file
func writeFragsToFastaFile(frags []*Frag, maxIDLength int, circularize bool, fastaFile *os.File) (err error) {
	truncID := func(s string) string {
		if len(s) < maxIDLength {
			return s
		} else {
			return s[:maxIDLength]
		}
	}

	// create a multimap of fragments indexed by truncated ID
	fragsByTruncatedIDs := make(map[string][]*Frag)
	for _, f := range frags {
		fragID := truncID(f.ID)
		fragsWithFragID, fragIDFound := fragsByTruncatedIDs[fragID]
		if fragIDFound {
			fragsByTruncatedIDs[fragID] = append(fragsWithFragID, f)
		} else {
			fragsByTruncatedIDs[fragID] = []*Frag{f}
		}
	}

	// convert an int to an Excel like column a ... z, aa .. az, ba .. bz
	// just to be safe
	base10ToBase26 := func(i int) string {
		var base26Val string = ""
		for currVal := i; ; {
			if currVal >= 26 {
				mod := currVal % 26
				currVal = currVal/26 - 1
				base26Val = fmt.Sprintf("%c", 'a'+rune(mod)) + base26Val
			} else {
				return fmt.Sprintf("%c", 'a'+rune(currVal)) + base26Val
			}
		}
	}

	for fragID, fragsWithFragID := range fragsByTruncatedIDs {
		if len(fragsWithFragID) == 1 {
			// no duplicates
			f := fragsWithFragID[0]
			rlog.Debugf("Write %s", f.ID)
			if ferr := writeSeqToFastaFile(fragID, f.Seq, circularize, fastaFile); ferr != nil {
				rlog.Errorf("Error writing fragment %s\n", f.ID)
				err = multierr.Append(err, ferr)
			}
		} else {
			// handle duplicates
			rlog.Infof("%d blast DB fragment ID duplicates found for %s", len(fragsWithFragID), fragID)
			for i, f := range fragsWithFragID {
				fragIDPrefix := fragIDComponents(f.ID)[0]
				fragIDSuffix := f.ID[len(fragIDPrefix):]
				newFragID := truncID(fmt.Sprintf("%s%s%s", fragIDPrefix, base10ToBase26(i), fragIDSuffix))

				if ferr := writeSeqToFastaFile(newFragID, f.Seq, circularize, fastaFile); ferr != nil {
					rlog.Errorf("Error writing fragment %s\n", f.ID)
					err = multierr.Append(err, ferr)
				}
			}
		}
	}

	return err
}

func writeSeqToFastaFile(id, seq string, circular bool, fastaFile *os.File) (err error) {
	var outputSeq, circularAttr string
	if circular {
		firstHalf := seq[:len(seq)/2]
		secondHalf := seq[len(seq)/2:]
		if firstHalf == secondHalf {
			// the current sequence is already circularized
			outputSeq = seq
		} else {
			outputSeq = seq + seq
		}
		circularAttr = "circular"
	} else {
		outputSeq = seq
		circularAttr = ""
	}
	_, err = fastaFile.WriteString(fmt.Sprintf(">%s %s\n%s\n", id, circularAttr, outputSeq))
	return err
}

// writeGenbank writes a slice of fragments/features to a genbank output file.
func writeGenbank(filename, name, seq string, frags []*Frag, feats []match) {
	// header row
	d := time.Now().Local()
	h1 := fmt.Sprintf("LOCUS       %s", name)
	h2 := fmt.Sprintf("%d bp DNA      circular      %s\n", len(seq), strings.ToUpper(d.Format("02-Jan-2006")))
	space := strings.Repeat(" ", 81-len(h1+h2))
	header := h1 + space + h2

	// feature rows
	var fsb strings.Builder
	fsb.WriteString("DEFINITION  .\nACCESSION   .\nFEATURES             Location/Qualifiers\n")
	for _, m := range feats {
		cS := ""
		cE := ""
		if m.isRevCompMatch() {
			cS = "complement("
			cE = ")"
		}

		s := (m.queryStart + 1) % len(seq)
		e := (m.queryEnd + 1) % len(seq)

		if s == 0 {
			s = len(seq)
		}
		if e == 0 {
			e = len(seq)
		}

		fsb.WriteString(
			fmt.Sprintf("     misc_feature    %s%d..%d%s\n", cS, s, e, cE) +
				fmt.Sprintf("                     /label=\"%s\"\n", m.entry),
		)
	}

	// origin row
	var ori strings.Builder
	ori.WriteString("ORIGIN\n")
	for i := 0; i < len(seq); i += 60 {
		n := strconv.Itoa(i + 1)
		ori.WriteString(strings.Repeat(" ", 9-len(n)) + n)
		for s := i; s < i+60 && s < len(seq); s += 10 {
			e := s + 10
			if e > len(seq) {
				e = len(seq)
			}
			ori.WriteString(fmt.Sprintf(" %s", seq[s:e]))
		}
		ori.WriteString("\n")
	}
	ori.WriteString("//\n")

	gb := strings.Join([]string{header, fsb.String(), ori.String()}, "")
	err := os.WriteFile(filename, []byte(gb), 0644)
	if err != nil {
		rlog.Fatal(err)
	}
}
