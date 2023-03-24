package repp

import (
	"encoding/csv"
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"
)

type oligo struct {
	id        string
	seq       string
	markAsNew bool
}

func (o oligo) isEmpty() bool {
	return o.seq == ""
}

func (o oligo) hasID() bool {
	return o.id != ""
}

func (o oligo) getIDOrNA(markNewOligo bool) string {
	if o.hasID() {
		if markNewOligo && o.markAsNew {
			return "*" + o.id
		} else {
			return o.id
		}
	} else {
		return "N/A"
	}
}

func (o *oligo) assignNewOligoID(id string) {
	o.id = id
	o.markAsNew = true
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
	idBase1, index1 := extractOligoIDComps(a[i].id)
	idBase2, index2 := extractOligoIDComps(a[j].id)

	idBaseCmpRes := strings.Compare(idBase1, idBase2)
	if idBaseCmpRes == 0 {
		return index1 < index2
	} else {
		return idBaseCmpRes < 0
	}
}

func (a sortedOligosByID) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

type oligosDB struct {
	indexedOligos map[string]oligo
	oligoIDBase   string
	nextOligoID   uint
}

// check if the provided sequence exists in the database
// if it exists it returns a full oligo (that has both the ID and the sequence set)
// otherwise the oligo has only the sequence filled in
// if the provided sequence is empty it returns an empty oligo (no ID and no sequence set)
func (oligos oligosDB) findOligo(seq string) oligo {
	if seq == "" {
		return oligo{}
	}
	o, found := oligos.indexedOligos[strings.ToUpper(seq)]
	if found {
		return o
	} else {
		return oligo{seq: seq}
	}
}

func (oligos oligosDB) getNewOligoID(newSeqIndex int) string {
	return fmt.Sprintf("%s%d", oligos.oligoIDBase, oligos.nextOligoID+uint(newSeqIndex))
}

func newOligos() *oligosDB {
	var o = &oligosDB{}
	o.indexedOligos = make(map[string]oligo)
	o.oligoIDBase = "oS" // default ID base
	o.nextOligoID = 1
	return o
}

func readOligos(manifest string) (oligos *oligosDB) {
	oligos = newOligos()

	if manifest == "" {
		return
	}

	f, err := os.Open(manifest)
	if err != nil {
		rlog.Warnf("Error opening oligos manifest %s: %v", manifest, err)
		return
	}
	defer f.Close()

	manifestReader := csv.NewReader(f)

	if err = readOligosFromCSV(manifestReader, oligos); err != nil {
		rlog.Warnf("Error parsing oligos manifest %s: %v", manifest, err)
	}

	return
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
			oligos.oligoIDBase = currentOligoIDBase
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
			id:  oligoIdField,
			seq: oligoSeqField, // put the original sequence field here as read from the file
		}
		oligos.indexedOligos[strings.ToUpper(oligoSequence)] = oligo
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
