package repp

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"os/exec"
	"strconv"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
	"go.uber.org/multierr"
)

// primer3 is a utility struct for executing primer3 to create primers on a fragment
type primer3 struct {
	// the target sequence
	seq string

	// input file
	in *os.File

	// output file
	out *os.File

	// path to primer3 executable
	primer3Exec string

	// path to primer3 config folder (with trailing separator)
	// this is the path to primer thermodynamic parameters
	primer3ConfDir string

	// configuaration
	config *config.Config
}

// newPrimer3 creates a primer3 struct from a fragment
func newPrimer3(seq string, conf *config.Config) primer3 {
	// try to get primer3 executable from env if set
	primer3Exec := getExecutable("PRIMER3_HOME", "bin", "primer3_core")
	return primer3{
		seq:            strings.ToUpper(seq),
		primer3Exec:    primer3Exec,
		primer3ConfDir: conf.GetPrimer3ConfigDir(),
		config:         conf,
	}
}

// input makes a primer3 input settings file and writes it to the filesystem
//
// the primers on this Frag should account for creating homology
// against the last Frag and the next Frag if there isn't enough
// existing homology to begin with (the two nodes should share ~50/50)
//
// returning the number of bp that have to be artifically added to the left and right primers
func (p *primer3) input(f, prev, next *Frag) (addLeft, addRight int, err error) {
	in, inErr := os.CreateTemp("", "primer3-in-*")
	out, outErr := os.CreateTemp("", "primer3-out-*")

	if inErr != nil || outErr != nil {
		return 0, 0, multierr.Append(inErr, outErr)
	}
	p.in = in
	p.out = out

	// adjust the Frag's start and end index in the event that there's too much homology
	// with the neighboring fragment
	p.shrink(prev, f, next)

	// calc the bps to add on the left and right side of this Frag
	addLeft = p.bpToAdd(prev, f)
	addRight = p.bpToAdd(f, next)

	start := f.start
	length := f.end - start + 1

	// check whether we have wiggle room on the left or right hand sides to move the
	// primers inward (let primer3 pick better primers)
	//
	// also adjust start and length in case there's TOO large an overhang and we need
	// to trim it in one direction or the other
	leftBuffer := p.buffer(prev.distTo(f))
	rightBuffer := p.buffer(f.distTo(next))

	if length-leftBuffer-rightBuffer < p.config.PcrMinFragLength {
		leftBuffer = 0
		rightBuffer = 0
	}

	// create the settings map from all instructions
	settings := p.settings(
		f.ID,
		start,
		length,
		leftBuffer,
		rightBuffer,
	)
	// write the settings to a buffer
	var fileBuffer bytes.Buffer
	for key, val := range settings {
		fmt.Fprintf(&fileBuffer, "%s=%s\n", key, val)
	}
	fileBuffer.WriteString("=") // required at file's end
	// then write them to the file
	if _, err = p.in.Write(fileBuffer.Bytes()); err != nil {
		return 0, 0, fmt.Errorf("failed to write primer3 input file %v: ", err)
	}
	return
}

// shrink adjusts the start and end of a Frag in the scenario where
// it excessively overlaps a neighboring fragment. For example, if there's
// 700bp of overlap, this will trim it back so we just PCR a subselection of
// the Frag and keep the overlap beneath the upper limit.
// Only the end if shrunk. Only shifting from right side.
// Otherwise, two neighboring fragments will both shrink and there won't be an overlap
func (p *primer3) shrink(last, f, next *Frag) *Frag {
	var shiftInLeft int
	var shiftInRight int

	if distRight := f.distTo(next); distRight < -p.config.FragmentsMaxHomology {
		// there's too much homology on the right side, we should move the Frag's end inward
		shiftInRight = (-distRight) - p.config.FragmentsMaxHomology
	}

	// make sure the fragment doesn't become less than the minimum length
	canShrink := (f.end-shiftInRight)-(f.start+shiftInLeft) > p.config.PcrMinFragLength &&
		len(f.Seq)-shiftInRight > shiftInLeft
	if canShrink {
		f.start += shiftInLeft
		f.end -= shiftInRight
		f.Seq = f.Seq[shiftInLeft : len(f.Seq)-shiftInRight]
	}

	return f
}

// bpToAdd returns the number of bp to add the end of the left Frag to create a junction
// with the right Frag
func (p *primer3) bpToAdd(left, right *Frag) int {
	if !left.couldOverlapViaPCR(right) {
		return 0 // we're going to synthesize there, don't add bp via PCR
	}

	if left.overlapsViaHomology(right) {
		return 0 // there is already enough overlap via PCR
	}

	bpDist := left.distTo(right) + 1 // if there's a gap
	if bpDist < 0 {
		bpDist = 0
	}

	// this Frag will add half the homology to the last fragment
	// eg: 5 bp distance leads to 2.5bp + ~10bp additonal
	// eg: -10bp distance leads to ~0 bp additional:
	// 		other Frag is responsible for all of it
	b := math.Ceil(float64(p.config.FragmentsMinHomology) / 2)

	return bpDist + int(b)
}

// buffer takes the dist from a one fragment to another and
// returns the length of the "buffer" in which the primers can be optimized (let primer3 pick)
//
// dist is positive if there's a gap between the start/end of a fragment and the start/end of
// the other and negative if they overlap
func (p *primer3) buffer(dist int) (buffer int) {
	if dist > 2*p.config.PcrPrimerMaxEmbedLength-p.config.FragmentsMinHomology {
		// we'll synthesize because the gap is so large, add 100bp of buffer
		return p.config.PcrBufferLength
	}

	if dist < -p.config.FragmentsMinHomology {
		// there's enough additonal overlap that we can move this FWD primer inwards
		// but only enough to ensure that there's still minHomology bp overlap
		// and only enough so we leave the neighbor space for primer optimization too
		return (-dist - p.config.FragmentsMinHomology) / 2
	}

	return 0
}

// settingsMap returns a new settings map for the primer3 config files
// can either use pick_cloning_primers mode, if the start and end primers' locations
// are fixed, or pick_primer_list mode if we're letting the primers shift and allowing
// primer3 to pick the best ones. One side may be free to move and the other not
func (p *primer3) settings(seqID string, start, length, leftBuffer, rightBuffer int) map[string]string {
	var strictPrimerSelection string
	if p.config.PcrPrimerUseStrictConstraints {
		strictPrimerSelection = "0"
	} else {
		strictPrimerSelection = "1"
	}

	// see primer3 manual or /vendor/primer3-2.4.0/settings_files/p3_th_settings.txt
	settings := map[string]string{
		"SEQUENCE_ID":                          seqID,
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH": p.primer3ConfDir,
		"PRIMER_NUM_RETURN":                    "1",
		"PRIMER_PICK_ANYWAY":                   strictPrimerSelection,
		"SEQUENCE_TEMPLATE":                    p.seq + p.seq,                             // TODO
		"PRIMER_MIN_SIZE":                      strconv.Itoa(p.config.PcrPrimerMinLength), // default 18
		"PRIMER_OPT_SIZE":                      strconv.Itoa(p.config.PcrPrimerOptimumLength),
		"PRIMER_MAX_SIZE":                      strconv.Itoa(p.config.PcrPrimerMaxLength),
		"PRIMER_EXPLAIN_FLAG":                  "1",
		"PRIMER_MIN_TM":                        fmt.Sprintf("%f", p.config.PcrPrimerMinTm),          // defaults to 57.0
		"PRIMER_MAX_TM":                        fmt.Sprintf("%f", p.config.PcrPrimerMaxTm),          // defaults to 63.0
		"PRIMER_MAX_HAIRPIN_TH":                fmt.Sprintf("%f", p.config.FragmentsMaxHairpinMelt), // defaults to 47.0
	}
	if p.config.PcrMaxHomopolymerLength > 0 {
		settings["PRIMER_MAX_POLY_X"] = fmt.Sprintf("%d", p.config.PcrMaxHomopolymerLength) // defaults to 5
	}
	if p.config.PcrPairMaxBindingScore > 0 {
		settings["PRIMER_PAIR_MAX_COMPL_ANY"] = fmt.Sprintf("%.2f", p.config.PcrPairMaxBindingScore) // defaults to 8.00
	}
	// if there is room to optimize, we let primer3 pick the best primers available
	// with a range on either side of the fragment's start
	// https://primer3.org/manual.html#SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
	if leftBuffer > 0 || rightBuffer > 0 {
		// V-ls  v-le               v-rs   v-re
		// ---------------------------------
		leftEnd := start + leftBuffer + p.config.PcrPrimerMaxLength
		rightStart := start + length - rightBuffer - p.config.PcrPrimerMaxLength
		excludeLength := rightStart - leftEnd

		if excludeLength >= 0 && excludeLength > p.config.PcrPrimerMaxLength {
			settings["PRIMER_TASK"] = "generic"
			settings["PRIMER_PICK_LEFT_PRIMER"] = "1"
			settings["PRIMER_PICK_INTERNAL_OLIGO"] = "0"
			settings["PRIMER_PICK_RIGHT_PRIMER"] = "1"

			// ugly undoing of the above in case only one side has buffer
			if leftBuffer == 0 {
				settings["SEQUENCE_FORCE_LEFT_START"] = strconv.Itoa(start)
				settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = fmt.Sprintf(",,%d,%d ;", rightStart, rightBuffer+p.config.PcrPrimerMaxLength)
			}
			if rightBuffer == 0 {
				settings["SEQUENCE_FORCE_RIGHT_START"] = strconv.Itoa(start + length)
				settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = fmt.Sprintf("%d,%d,, ;", start, leftBuffer+p.config.PcrPrimerMaxLength)
			}

			if settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] == "" {
				settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = fmt.Sprintf("%d,%d,%d,%d ;", start, leftBuffer+p.config.PcrPrimerMaxLength, rightStart, rightBuffer+p.config.PcrPrimerMaxLength)
			}
			settings["PRIMER_PRODUCT_SIZE_RANGE"] = fmt.Sprintf("%d-%d", excludeLength, length)
		}
	}

	// otherwise force the start and end of the PCR range
	if settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] == "" {
		settings["PRIMER_TASK"] = "pick_cloning_primers"
		settings["SEQUENCE_INCLUDED_REGION"] = fmt.Sprintf("%d,%d", start, length)
		settings["PRIMER_PRODUCT_SIZE_RANGE"] = fmt.Sprintf("%d-%d", length, length)
	}

	return settings
}

// run the primer3 executable against the input file
func (p *primer3) run() (err error) {
	p3Cmd := exec.Command(
		p.primer3Exec,
		p.in.Name(),
		"-output", p.out.Name(),
		"-strict_tags",
	)

	// execute primer3 and wait on it to finish
	if output, err := p3Cmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to execute primer3 on input file %s: %s: %v", p.in.Name(), string(output), err)
	}

	return
}

// parse the output into primers. add to fragment
//
// target is the target sequence we're building for. We need it to modulo the primer ranges
func (p *primer3) parse(target string) (primers []Primer, err error) {
	fileBytes, err := os.ReadFile(p.out.Name())
	if err != nil {
		return
	}
	file := string(fileBytes)

	// read in results into map, they're all 1:1
	results := make(map[string]string)
	for _, line := range strings.Split(file, "\n") {
		keyVal := strings.Split(line, "=")
		if len(keyVal) > 1 {
			results[strings.TrimSpace(keyVal[0])] = strings.TrimSpace(keyVal[1])
		}
	}

	if p3Warnings := results["PRIMER_WARNING"]; p3Warnings != "" {
		err = fmt.Errorf("warnings executing primer3: %s", p3Warnings)
		return
	}

	if p3Error := results["PRIMER_ERROR"]; p3Error != "" {
		err = fmt.Errorf("failed to execute primer3 against %s: %s", file, p3Error)
		return
	}

	if numPairs := results["PRIMER_PAIR_NUM_RETURNED"]; numPairs == "0" {
		err = fmt.Errorf("failed to create primers using: \n%s", file)
		return
	}

	// read in a single primer from the output string file
	// side is either "LEFT" or "RIGHT"
	parsePrimer := func(side string, index int) Primer {
		seq := results[fmt.Sprintf("PRIMER_%s_%d_SEQUENCE", side, index)]
		tm := results[fmt.Sprintf("PRIMER_%s_%d_TM", side, index)]
		gc := results[fmt.Sprintf("PRIMER_%s_%d_GC_PERCENT", side, index)]
		penalty := results[fmt.Sprintf("PRIMER_%s_%d_PENALTY", side, index)]
		pairPenalty := results[fmt.Sprintf("PRIMER_PAIR_%d_PENALTY", index)]
		notes := results[fmt.Sprintf("PRIMER_%s_%d_PROBLEMS", side, index)]

		tmValue, _ := strconv.ParseFloat(tm, 64)
		gcValue, _ := strconv.ParseFloat(gc, 64)
		penaltyValue, _ := strconv.ParseFloat(penalty, 64)
		pairValue, _ := strconv.ParseFloat(pairPenalty, 64)

		primerRange := results[fmt.Sprintf("PRIMER_%s_0", side)]
		primerStart, _ := strconv.Atoi(strings.Split(primerRange, ",")[0])
		primerEnd := primerStart + len(seq)
		if side == "RIGHT" {
			primerStart -= len(seq)
			primerEnd = primerStart + len(seq)
		}

		return Primer{
			Seq:           seq,
			Strand:        side == "LEFT",
			Tm:            tmValue,
			GC:            gcValue,
			Penalty:       penaltyValue,
			PairPenalty:   pairValue,
			PrimingRegion: seq,
			Range: ranged{
				start: primerStart,
				end:   primerEnd,
			},
			Notes: notes,
		}
	}
	primers = []Primer{
		parsePrimer("LEFT", 0),
		parsePrimer("RIGHT", 0),
	}
	return
}

func (p *primer3) close() (err error) {
	if os.Getenv("DEBUG_REPP") == "TRUE" {
		// keep the temporary files
		rlog.Infof("Primer3 input/output: %s, %s", p.in.Name(), p.out.Name())
		return
	}
	// remove temporary input and output
	if p.in != nil {
		err = multierr.Append(err, os.Remove(p.in.Name()))
	}
	if p.out != nil {
		err = multierr.Append(err, os.Remove(p.out.Name()))
	}
	return
}

// hairpin finds the melting temperature of a hairpin in a sequence
// returns 0 if there is none
func hairpin(seq string, conf *config.Config) (melt float64) {
	// if it's longer than 60bp (max for ntthal) find the max between
	// the start and end of the sequence
	if len(seq) > 60 {
		startHairpin := hairpin(seq[:60], conf)
		endHairpin := hairpin(seq[len(seq)-60:], conf)

		if startHairpin > endHairpin {
			return startHairpin
		}
		return endHairpin
	}

	// see nnthal (no parameters) help. within primer3 distribution
	ntthalCmd := exec.Command(
		getExecutable("PRIMER3_HOME", "bin", "ntthal"),
		"-a", "HAIRPIN",
		"-r",       // temperature only
		"-t", "50", // gibson assembly is at 50 degrees
		"-s1", seq,
		"-path", conf.GetPrimer3ConfigDir(),
	)

	ntthalOut, err := ntthalCmd.CombinedOutput()
	if err != nil {
		stderr.Printf("failed to execute ntthal: -s1 %s -path %s", seq, conf.GetPrimer3ConfigDir())
		rlog.Fatal(err)
	}

	ntthalOutString := string(ntthalOut)
	temp, err := strconv.ParseFloat(strings.TrimSpace(ntthalOutString), 64)
	if err != nil {
		stderr.Printf("failed to parse ntthal: -s1 %s -path %s", seq, conf.GetPrimer3ConfigDir())
		rlog.Fatal(err)
	}

	return temp
}

// reverseComplement returns the reverse complement of a sequence
func reverseComplement(seq string) string {
	seq = strings.ToUpper(seq)

	revCompMap := map[rune]byte{
		'A': 'T',
		'T': 'A',
		'G': 'C',
		'C': 'G',
		'^': '_',
		'_': '^',
	}

	var revCompBuffer bytes.Buffer
	for _, c := range seq {
		revCompBuffer.WriteByte(revCompMap[c])
	}

	revCompBytes := revCompBuffer.Bytes()
	for i := 0; i < len(revCompBytes)/2; i++ {
		j := len(revCompBytes) - i - 1
		revCompBytes[i], revCompBytes[j] = revCompBytes[j], revCompBytes[i]
	}

	return string(revCompBytes)
}
