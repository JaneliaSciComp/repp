package repp

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"go.uber.org/multierr"
)

var (
	// stderr is for logging to Stderr (without an annoying timestamp)
	stderr = log.New(os.Stderr, "", 0)
)

type AssemblyParams interface {
	GetIn() string
	SetIn(in string)

	GetOut() string
	SetOut(out string)

	GetOutputFormat() string
	SetOutputFormat(f string)

	GetFilters() []string
	SetFilters(fs []string)

	GetIdentity() int
	SetIdentity(i int)

	GetBackboneName() string
	SetBackboneName(bn string)

	getDBs() ([]DB, error)
	SetDbNames(dbNames []string)

	getEnzymes() ([]enzyme, error)
	SetEnzymeNames(enzymeNames []string)
}

// assemblyParamsImpl contains assembly input parameters.
type assemblyParamsImpl struct {
	// the name of the file to read the input from
	in string

	// the name of the file to write the output to
	out string

	// output format (JSON, CSV)
	outFormat string

	// a list of dbs to run BLAST against (their names' on the filesystem)
	dbNames []string

	// the backbone (optional) to insert the pieces into
	backboneName string

	// list of enzimes
	enzymeNames []string

	// slice of strings to weed out fragments from BLAST matches
	filters []string

	// percentage identity for finding building fragments in BLAST databases
	identity int
}

func MkAssemblyParams() AssemblyParams {
	return &assemblyParamsImpl{}
}

func (ap assemblyParamsImpl) GetIn() string {
	return ap.in
}

func (ap *assemblyParamsImpl) SetIn(in string) {
	ap.in = in
}

func (ap assemblyParamsImpl) GetOut() string {
	return ap.out
}

func (ap *assemblyParamsImpl) SetOut(out string) {
	ap.out = out
}

func (ap assemblyParamsImpl) GetOutputFormat() string {
	return ap.outFormat
}

func (ap *assemblyParamsImpl) SetOutputFormat(f string) {
	ap.outFormat = f
}

func (ap assemblyParamsImpl) GetFilters() []string {
	return ap.filters
}

func (ap *assemblyParamsImpl) SetFilters(filters []string) {
	ap.filters = filters
}

func (ap assemblyParamsImpl) GetIdentity() int {
	return ap.identity
}

func (ap *assemblyParamsImpl) SetIdentity(identity int) {
	ap.identity = identity
}

func (ap assemblyParamsImpl) GetBackboneName() string {
	return ap.backboneName
}

func (ap *assemblyParamsImpl) SetBackboneName(backboneName string) {
	ap.backboneName = backboneName
}

func (ap assemblyParamsImpl) getDBs() (dbs []DB, err error) {
	return getRegisteredDBs(ap.dbNames)
}

func (ap *assemblyParamsImpl) SetDbNames(dbNames []string) {
	ap.dbNames = dbNames
}

func (ap assemblyParamsImpl) getEnzymes() (enzymes []enzyme, err error) {
	return getValidEnzymes(ap.enzymeNames)
}

func (ap *assemblyParamsImpl) SetEnzymeNames(enzymeNames []string) {
	ap.enzymeNames = enzymeNames
}

func prepareBackbone(
	bbName string,
	enzymes []enzyme,
	dbs []DB) (f *Frag, backbone *Backbone, err error) {

	if bbName == "" {
		// if no backbone was specified, return an empty Frag
		return &Frag{}, &Backbone{}, nil
	}

	// confirm that the backbone exists in one of the dbs (or local fs) gather it as a Frag if it does
	bbFrag, err := queryDatabases(bbName, dbs)
	if err != nil {
		return &Frag{}, &Backbone{}, err
	}

	// try to digest the backbone with the enzyme
	if len(enzymes) == 0 {
		return &Frag{},
			&Backbone{},
			fmt.Errorf("backbone passed, %s, without an enzyme to digest it", bbName)
	}

	if f, backbone, err = digest(bbFrag, enzymes); err != nil {
		return &Frag{}, &Backbone{}, err
	}

	return
}

// read a dir of FASTA or Genbank files to a slice of fragments
func multiFileRead(fs []string) (fragments []*Frag, err error) {
	for _, f := range fs {
		fFrags, ferr := read(f, false)
		if ferr != nil {
			rlog.Warnf("Error reading sequence from %s\n", f)
			err = multierr.Append(err, ferr)
		} else {
			fragments = append(fragments, fFrags...)
		}
	}

	return
}

// read a FASTA or Genbank file (by its path on local FS) to a slice of Fragments.
func read(path string, feature bool) (fragments []*Frag, err error) {
	if !filepath.IsAbs(path) {
		path, err = filepath.Abs(path)
		if err != nil {
			return nil, fmt.Errorf("failed to create path to input file: %s", err)
		}
	}

	dat, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	file := string(dat)

	path = strings.ToLower(path)
	if strings.HasSuffix(path, "fa") ||
		strings.HasSuffix(path, "fasta") ||
		file[0] == '>' {
		return readFasta(path, file)
	}

	if strings.HasSuffix(path, "gb") ||
		strings.HasSuffix(path, "gbk") ||
		strings.HasSuffix(path, "genbank") {
		return readGenbank(path, file, feature)
	}

	return nil, fmt.Errorf("failed to parse %s: unrecognized file type", path)
}

// readFasta parses the multifasta file to fragments.
func readFasta(path, contents string) (frags []*Frag, err error) {
	// split by newlines
	lines := strings.Split(contents, "\n")

	// read in the frags
	var headerIndices []int
	var ids []string
	var fragTypes []fragType
	for i, line := range lines {
		if strings.HasPrefix(line, ">") {
			headerIndices = append(headerIndices, i)
			ids = append(ids, strings.TrimSpace(line[1:]))
			if strings.Contains(line, "circular") {
				fragTypes = append(fragTypes, circular)
			} else {
				fragTypes = append(fragTypes, linear)
			}
		}
	}

	// create a regex for cleaning the sequence
	var unwantedChars = regexp.MustCompile(`(?im)[^atgc]|\W`)

	// accumulate the sequences from between the headers
	var seqs []string
	for i, headerIndex := range headerIndices {
		nextLine := len(lines)
		if i < len(headerIndices)-1 {
			nextLine = headerIndices[i+1]
		}
		seqLines := lines[headerIndex+1 : nextLine]
		seqJoined := strings.Join(seqLines, "")
		seq := unwantedChars.ReplaceAllString(seqJoined, "")
		seq = strings.ToUpper(seq)
		seqs = append(seqs, seq)
	}

	// build and return the new frags
	for i, id := range ids {
		frags = append(frags, &Frag{
			ID:       id,
			Seq:      seqs[i],
			fragType: fragTypes[i],
		})
	}

	// opened and parsed file but found nothing
	if len(frags) < 1 {
		return frags, fmt.Errorf("failed to parse fragment(s) from %s", path)
	}

	return
}

// readGenbank parses a genbank file to fragments. Returns either fragments or parseFeatures,
// depending on the parseFeatures parameter.
func readGenbank(path, contents string, parseFeatures bool) (fragments []*Frag, err error) {
	genbankSplit := strings.Split(contents, "ORIGIN")

	if len(genbankSplit) != 2 {
		return nil, fmt.Errorf("failed to parse %s: improperly formatted genbank file", path)
	}

	seq := strings.ToUpper(genbankSplit[1])
	nonBpRegex := regexp.MustCompile("[^ATGC]")
	cleanedSeq := nonBpRegex.ReplaceAllString(seq, "")

	if parseFeatures {
		// parse each feature to a fragment (misnomer)
		splitOnFeatures := strings.Split(genbankSplit[0], "FEATURES")

		if len(splitOnFeatures) < 2 {
			return nil, fmt.Errorf("failed to parse features from %s", path)
		}

		featureSplitRegex := regexp.MustCompile(`\w+\s+\w+`)
		featureStrings := featureSplitRegex.Split(splitOnFeatures[1], -1)

		features := []*Frag{}
		for featureIndex, feature := range featureStrings {
			rangeRegex := regexp.MustCompile(`(\d*)\.\.(\d*)`)
			rangeIndexes := rangeRegex.FindStringSubmatch(feature)

			if len(rangeIndexes) < 3 {
				continue
			}

			start, err := strconv.Atoi(rangeIndexes[1])
			if err != nil {
				return nil, err
			}

			end, err := strconv.Atoi(rangeIndexes[2])
			if err != nil {
				return nil, err
			}
			featureSeq := cleanedSeq[start-1 : end] // make 0-indexed
			featureSeq = strings.ToUpper(featureSeq)

			labelRegex := regexp.MustCompile(`\/label=(.*)`)
			labelMatch := labelRegex.FindStringSubmatch(feature)
			label := ""
			if len(labelMatch) > 1 {
				label = labelMatch[1]
			} else {
				label = strconv.Itoa(featureIndex)
			}

			features = append(features, &Frag{
				ID:  label,
				Seq: featureSeq,
			})
		}

		return features, nil
	}

	// parse just the file's sequence
	idRegex := regexp.MustCompile(`LOCUS[ \t]*([^ \t]*)`)
	idMatches := idRegex.FindStringSubmatch(genbankSplit[0])

	var id string
	if len(idMatches) == 0 {
		return nil, fmt.Errorf("failed to parse locus from %s", path)
	} else if len(idMatches) > 1 {
		id = idMatches[1]
	} else {
		// use filename otherwise if the ID is just LOCUS
		// and if other files have that there will be bad ids
		id = filepath.Base(path)
	}

	return []*Frag{
		{
			ID:  id,
			Seq: cleanedSeq,
		},
	}, nil
}
