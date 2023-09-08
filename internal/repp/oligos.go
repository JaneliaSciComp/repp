package repp

import (
	"encoding/csv"
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"

	"go.uber.org/multierr"
)

type oligo struct {
	id    string
	seq   string
	isNew bool
	synth bool
}

func (o oligo) isEmpty() bool {
	return o.seq == ""
}

func (o oligo) hasID() bool {
	return o.id != ""
}

func (o oligo) getIDOrDefault(markID bool, defaultValue string) string {
	if o.hasID() {
		if markID {
			return "*" + o.id
		} else {
			return o.id
		}
	} else {
		return defaultValue
	}
}

func (o *oligo) assignNewOligoID(id string) {
	o.id = id
	o.isNew = true
}

var (
	oligoIDPattern *regexp.Regexp = regexp.MustCompile(`(?P<Base>\D*)(?P<Index>\d*)$`)
	baseMatchPos   int            = oligoIDPattern.SubexpIndex("Base")
	indexMatchPos  int            = oligoIDPattern.SubexpIndex("Index")
)

type sortedOligosByID []oligo

func (a sortedOligosByID) Len() int {
	return len(a)
}

func (a sortedOligosByID) Less(i, j int) bool {
	if a[i].synth && !a[j].synth {
		return false
	} else if !a[i].synth && a[i].synth {
		return true
	} else if !a[i].synth && !a[j].synth {
		return compOligoIDs(a[i], a[j]) < 0
	} else {
		return strings.Compare(a[i].id, a[j].id) < 0
	}
}

func compOligoIDs(o1, o2 oligo) int {
	idBase1, index1 := extractOligoIDComps(o1.id)
	idBase2, index2 := extractOligoIDComps(o2.id)

	idBaseCmpRes := strings.Compare(idBase1, idBase2)
	if idBaseCmpRes == 0 {
		if index1 < index2 {
			return -1
		} else if index1 > index2 {
			return 1
		} else {
			return 0
		}
	} else {
		return idBaseCmpRes
	}
}

func (a sortedOligosByID) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

type oligosDB struct {
	indexedOligos     map[string]oligo
	oligoIDBasePrefix string
	nextOligoID       uint
	synthOligos       bool
}

func (oligos oligosDB) getNewOligoID(newSeqIndex int) string {
	// otherwise use the oligoIDBase as a prefix
	return fmt.Sprintf("%s%d", oligos.oligoIDBasePrefix, oligos.nextOligoID+uint(newSeqIndex))
}

func newOligosDB(defaultBasePrefix string, synthOligos bool) *oligosDB {
	var o = &oligosDB{}
	o.indexedOligos = make(map[string]oligo)
	o.oligoIDBasePrefix = defaultBasePrefix
	o.nextOligoID = 1
	o.synthOligos = synthOligos
	return o
}

func (oligos *oligosDB) addOligo(o oligo) {
	oligos.indexedOligos[strings.ToUpper(o.seq)] = o
}

// check if the provided sequence exists in the provided databases
// if it exists it returns a full oligo (that has both the ID and the sequence set)
// otherwise the oligo has only the sequence filled in
// if the provided sequence is empty it returns an empty oligo (no ID and no sequence set)
func searchOligoDBs(seq string, oligoDBs []*oligosDB) oligo {
	if seq == "" {
		return oligo{}
	}
	for _, oligoDB := range oligoDBs {
		o, found := oligoDB.indexedOligos[strings.ToUpper(seq)]
		if found {
			return o
		}
	}
	return oligo{seq: seq}
}

func readOligos(dbLocations []string, basePrefix string, synthOligos bool) (oligos *oligosDB) {
	oligos = newOligosDB(basePrefix, synthOligos)
	oligosFnames, collectFilesErr := CollectFiles(dbLocations)
	if collectFilesErr != nil {
		rlog.Warnf("Errors trying to collect oligo filenames from: %v", dbLocations)
		return
	}

	var allErrs error
	for _, fn := range oligosFnames {
		if fn == "" {
			continue
		}
		oligosReadErr := readOligosFromFile(fn, oligos)
		if oligosReadErr != nil {
			allErrs = multierr.Append(allErrs, oligosReadErr)
		}
	}

	if allErrs != nil {
		rlog.Warnf("Errors trying to read oligos: %v", allErrs)
	}
	return
}

func readOligosFromFile(oligosCSVFilename string, oligos *oligosDB) error {
	f, err := os.Open(oligosCSVFilename)
	if err != nil {
		rlog.Warnf("Error opening oligos manifest %s: %v", oligosCSVFilename, err)
		return err
	}
	defer f.Close()

	manifestReader := csv.NewReader(f)

	if err = readOligosFromCSV(manifestReader, oligos); err != nil {
		rlog.Warnf("Error parsing oligos manifest %s: %v", oligosCSVFilename, err)
	}

	return nil
}

func readOligosFromCSV(manifestReader *csv.Reader, oligos *oligosDB) error {

	manifestReader.Comment = '#'
	manifestReader.TrimLeadingSpace = true
	records, err := manifestReader.ReadAll()
	if err != nil {
		return err
	}

	for i, r := range records {
		if len(r) < 2 {
			// skip this row because it has too few items
			rlog.Warnf("Skip row %d:%v because it has too few colums\n", i+1, r)
			continue
		}
		oligoIdField := strings.TrimSpace(r[0])
		oligoSeqField := strings.TrimSpace(r[1])

		if strings.EqualFold(oligoSeqField, "sequence") {
			// this is header
			continue
		} else if oligoIdField == "" || oligoSeqField == "" {
			rlog.Warnf("Skip row %d:%v because ID and/or sequence field is empty\n", i+1, r)
			continue
		}

		currentOligoIDBase, currentOligoIndex := extractOligoIDComps(oligoIdField)
		if currentOligoIDBase != "" {
			oligos.oligoIDBasePrefix = currentOligoIDBase
		}
		// if an index is found - nextIndex becomes current index +1
		// otherwise increment previous value
		if currentOligoIndex == 0 {
			oligos.nextOligoID++
		} else {
			oligos.nextOligoID = currentOligoIndex + 1
		}

		// there are cases where the sequence field looks like "/5Biosg/GTGAAGTTCCCAAAGGTGCA"
		// so in this case I want to use only the nucleotide sequence for indexing
		sequenceStart := strings.LastIndex(oligoSeqField, "/")
		var oligoSequence string
		if sequenceStart == -1 {
			oligoSequence = oligoSeqField
		} else {
			oligoSequence = oligoSeqField[sequenceStart+1:]
		}
		oligo := oligo{
			id:    oligoIdField,
			seq:   oligoSequence, // put the original sequence field here as read from the file
			synth: oligos.synthOligos,
		}
		oligos.addOligo(oligo)
	}

	return nil
}

func extractOligoIDComps(oligoId string) (string, uint) {
	oligoIDMatch := oligoIDPattern.FindStringSubmatch(oligoId)
	if oligoIDMatch == nil {
		// if the ID does not match the expected pattern simply increment the ID
		return oligoId, 0
	} else {
		oligoIDBase := oligoIDMatch[baseMatchPos]
		oligoIndex := oligoIDMatch[indexMatchPos]
		// if an index is found - nextIndex becomes current index +1
		// otherwise increment previous nextIndex
		if oligoIndex != "" {
			oligoIndexVal, err := strconv.ParseUint(oligoIndex, 0, 32)
			if err != nil {
				return oligoIDBase, 0
			} else {
				return oligoIDBase, uint(oligoIndexVal)
			}
		} else {
			return oligoIDBase, 0
		}
	}
}
