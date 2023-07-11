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

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`
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
	oligos *oligosDB,
	insertSeqLength int,
	seconds float64,
	backbone *Backbone,
	conf *config.Config,
) (*Output, error) {
	out, err := prepareSolutionsOutput(
		targetName,
		targetSeq,
		assemblies,
		insertSeqLength,
		seconds,
		backbone,
		conf,
	)
	if err != nil {
		return nil, err
	}
	if format == "CSV" {
		err = writeCSV(filename, fragmentBase(filename), oligos, out)
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
	insertSeqLength int,
	seconds float64,
	backbone *Backbone,
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
		assemblyFragmentIDs := make(map[string]bool)
		gibson := false // whether it will be assembled via Gibson assembly
		hasPCR := false // whether there will be a batch PCR

		for _, f := range assembly {
			if f.fragType != linear && f.fragType != circular {
				gibson = true
			}

			if f.fragType == pcr {
				hasPCR = true
			}

			f.Type = f.fragType.String() // freeze fragment type

			// round to two decimal places
			if f.Cost, err = roundCost(f.cost(true)); err != nil {
				return nil, err
			}

			// if it's already in the assembly, don't count cost twice
			if _, contained := assemblyFragmentIDs[f.ID]; f.ID != "" && contained {
				if f.Cost, err = roundCost(f.cost(false)); err != nil {
					return nil, err // ignore repo procurement costs
				}
			} else {
				assemblyFragmentIDs[f.ID] = true
			}

			// accumulate assembly cost
			assemblyCost += f.Cost
		}

		if gibson {
			assemblyCost += conf.GibsonAssemblyCost + conf.GibsonAssemblyTimeCost
		}

		if hasPCR {
			assemblyCost += conf.PcrTimeCost
		}

		solutionCost, err := roundCost(assemblyCost)
		if err != nil {
			return nil, err
		}

		solutions = append(solutions, Solution{
			Count:     len(assembly),
			Cost:      solutionCost,
			Fragments: assembly,
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
func writeCSV(filename, fragmentIDBase string, oligos *oligosDB, out *Output) (err error) {

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
	err = strategyCSVWriter.Write([]string{
		"Frag ID",
		"Fwd Primer",
		"Rev Primer",
		"Template",
		"Size",
		"Match Pct",
	})
	if err != nil {
		return nil
	}
	// Write the reagents headers
	err = reagentsCSVWriter.Write([]string{
		"Reagent ID",
		"Seq",
	})
	for si, s := range out.Solutions {
		snumber := si + 1
		// Write the solution cost and the number of fragments
		if _, err = fmt.Fprintf(strategyFile,
			"# Solution %d\n# Fragments:%d, Cost: %f\n",
			snumber,
			s.Count, s.Cost); err != nil {
			return err
		}
		if _, err = fmt.Fprintf(reagentsFile, "# Solution %d\n", snumber); err != nil {
			return err
		}
		reagents := []oligo{}
		var newOligoSequenceIndex int = 0
		newOligos := newOligosDB()
		var allOligoDBs []*oligosDB = []*oligosDB{
			oligos,
			newOligos,
		}
		for fi, f := range s.Fragments {
			fnumber := fi + 1
			var fwdPrimerSeq, revPrimerSeq, synthSeq string

			fID := fmt.Sprintf("%s_%d_%s", fragmentIDBase, fnumber, fragTypeAsString(f.fragType))
			fwdPrimerSeq = f.getPrimerSeq(true)
			revPrimerSeq = f.getPrimerSeq(false)
			if fwdPrimerSeq == "" && revPrimerSeq == "" {
				synthSeq = f.Seq
			}

			fwdOligo := searchOligoDBs(fwdPrimerSeq, allOligoDBs)
			if !fwdOligo.isEmpty() {
				if !fwdOligo.hasID() {
					fwdOligo.assignNewOligoID(oligos.getNewOligoID(fragmentIDBase, newOligoSequenceIndex))
					newOligos.addOligo(fwdOligo)
					newOligoSequenceIndex++
				}
				reagents = append(reagents, fwdOligo)
			}
			revOligo := searchOligoDBs(revPrimerSeq, allOligoDBs)
			if !revOligo.isEmpty() {
				if !revOligo.hasID() {
					revOligo.assignNewOligoID(oligos.getNewOligoID(fragmentIDBase, newOligoSequenceIndex))
					newOligos.addOligo(revOligo)
					newOligoSequenceIndex++
				}
				reagents = append(reagents, revOligo)
			}
			var templateID string
			var matchRatio string
			if f.fragType == synthetic {
				templateID = "N/A"
				matchRatio = "N/A"
				synthReagent := oligo{
					id:    fID,
					seq:   synthSeq,
					synth: true,
				}
				reagents = append(reagents, synthReagent)
			} else {
				templateID = fragmentBase(f.ID)
				matchRatio = fmt.Sprintf("%d", int(f.matchRatio*100))
			}
			if err = strategyCSVWriter.Write([]string{
				fID,
				fwdOligo.getIDOrDefault(false, "N/A"), // fwd primer
				revOligo.getIDOrDefault(false, "N/A"), // rev primer
				templateID,                            // template
				strconv.Itoa(len(f.Seq)),
				matchRatio,
			}); err != nil {
				return nil
			}
		}
		strategyCSVWriter.Flush()
		sort.Sort(sortedOligosByID(reagents))
		for _, r := range reagents {
			reagentID := r.getIDOrDefault(!r.isNew && !r.synth, "N/A") // mark the ID if this reagent already existed in the original manifest
			err = writeReagent(reagentsCSVWriter, reagentID, r.seq)
			if err != nil {
				rlog.Errorf("Error writing reagent %s: %v", reagentID, err)
			}
		}
		reagentsCSVWriter.Flush()
	}

	return nil
}

func fragmentBase(filename string) string {
	var fragmentIDSeparators = " ,_-."

	splitFunc := func(c rune) bool {
		return strings.ContainsRune(fragmentIDSeparators, c)
	}

	baseNameFromFilename := strings.FieldsFunc(filepath.Base(filename), splitFunc)[0]

	if len(baseNameFromFilename) > 10 {
		// truncate if name is too long
		return baseNameFromFilename[:10]
	} else {
		return baseNameFromFilename
	}
}

func resultFilename(template, suffix string) string {
	ext := filepath.Ext(template)
	noExt := template[0 : len(template)-len(ext)]
	return noExt + "-" + suffix + ext
}

func writeReagent(csvWriter *csv.Writer, reagentID, reagentSeq string) (err error) {
	if reagentID != "" {
		err = csvWriter.Write([]string{
			reagentID,
			reagentSeq,
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
func writeFragsToFastaFile(frags []*Frag, fastaFile *os.File) (err error) {
	for _, f := range frags {
		rlog.Debugf("Write %s", f.ID)
		if _, ferr := fastaFile.WriteString(fmt.Sprintf(">%s\n%s\n", f.ID, f.Seq)); ferr != nil {
			rlog.Errorf("Error writing fragment %s\n", f.ID)
			err = multierr.Append(err, ferr)
		}
	}

	return
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
		if !m.forward {
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
