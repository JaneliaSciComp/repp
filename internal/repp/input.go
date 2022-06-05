package repp

import (
	"errors"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/spf13/cobra"
)

var (
	// stderr is for logging to Stderr (without an annoying timestamp)
	stderr = log.New(os.Stderr, "", 0)
)

// Flags contains parsed cobra Flags like "in", "out", "dbs", etc that are used by multiple commands.
type Flags struct {
	// the name of the file to write the input from
	in string

	// the name of the file to write the output to
	out string

	// a list of dbs to run BLAST against (their names' on the filesystem)
	dbs []DB

	// the backbone (optional) to insert the pieces into
	backbone *Frag

	// backbone meta (name, enzyme used to cut it, cut index)
	backboneMeta *Backbone

	// slice of strings to weed out fragments from BLAST matches
	filters []string

	// percentage identity for finding building fragments in BLAST databases
	identity int
}

// inputParser contains methods for parsing flags from the input &cobra.Command.
type inputParser struct{}

// NewFlags makes a new flags object manually. for testing.
func NewFlags(
	in, out, backbone, filter string,
	enzymes []string,
	dbs []DB,
) (*Flags, *config.Config) {
	c := config.New()

	p := inputParser{}
	parsedBB, bbMeta, err := p.parseBackbone(backbone, enzymes, dbs, c)
	if err != nil {
		rlog.Fatal(err)
	}

	if strings.Contains(in, ",") {
		in = p.parseFeatureInput(strings.Fields(in))
	}

	return &Flags{
		in:           in,
		out:          out,
		dbs:          dbs,
		backbone:     parsedBB,
		backboneMeta: bbMeta,
		filters:      p.getFilters(filter),
		identity:     98,
	}, c
}

// parseCmdFlags gathers the in path, out path, etc from a cobra cmd object
// returns Flags and a Config struct for repp.Plasmid or repp.Fragments.
func parseCmdFlags(cmd *cobra.Command, args []string, strict bool) (*Flags, *config.Config) {
	cmdName := strings.ToLower(cmd.Name())

	var err error
	fs := &Flags{} // parsed flags
	p := inputParser{}
	c := config.New()

	if fs.in, err = cmd.Flags().GetString("in"); fs.in == "" || err != nil {
		if cmdName == "features" {
			fs.in = p.parseFeatureInput(args)
		} else if cmdName == "sequence" && len(args) > 0 {
			fs.in = "input.fa"
			if err = ioutil.WriteFile(fs.in, []byte(fmt.Sprintf(">target_sequence\n%s", args[0])), 0644); err != nil {
				rlog.Fatal(err)
			}
		} else if fs.in, err = p.guessInput(); strict && err != nil {
			// check whether an input fail was specified
			if helperr := cmd.Help(); helperr != nil {
				rlog.Fatal(helperr)
			}
			rlog.Fatal(err)
		}
	}

	if fs.out, err = cmd.Flags().GetString("out"); strict && (fs.out == "" || err != nil) {
		fs.out = p.guessOutput(fs.in) // guess at an output name

		if fs.out == "" {
			if helperr := cmd.Help(); helperr != nil {
				rlog.Fatal(helperr)
			}
			rlog.Fatal("no output path")
		}
	}

	filters, err := cmd.Flags().GetString("exclude")
	if strict && err != nil && cmdName != "fragments" {
		if helperr := cmd.Help(); helperr != nil {
			rlog.Fatal(helperr)
		}
		rlog.Fatal("failed to parse filters: %v", err)
	}
	// try to split the filter fields into a list
	fs.filters = p.getFilters(filters)

	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = 100 // might be something other than `repp plasmid`
	}
	// set identity for blastn searching
	fs.identity = identity

	// read in the BLAST DB paths
	dbString, err := cmd.Flags().GetString("dbs")
	if err != nil {
		rlog.Fatal("failed to get dbs flag: %v", err)
	}
	m, err := newManifest()
	if err != nil {
		rlog.Fatalf("failed to get DB manifest: %v", err)
	}
	if fs.dbs, err = p.parseDBs(m, dbString); err != nil || len(fs.dbs) == 0 {
		rlog.Fatal("failed to find any fragment databases: %v", err)
	}

	// check if user asked for a specific backbone, confirm it exists in one of the dbs
	backbone, _ := cmd.Flags().GetString("backbone")

	// check if they also specified an enzyme
	enzymeList, _ := cmd.Flags().GetString("enzymeList")
	enzymes := p.parseCommaList(enzymeList)

	// try to digest the backbone with the enzyme
	fs.backbone, fs.backboneMeta, err = p.parseBackbone(backbone, enzymes, fs.dbs, c)
	if strict && err != nil {
		rlog.Fatal(err)
	}

	return fs, c
}

// guessInput returns the first fasta file in the current directory. Is used
// if the user hasn't specified an input file.
func (p *inputParser) guessInput() (in string, err error) {
	dir, _ := filepath.Abs(".")
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return
	}

	for _, file := range files {
		if file.IsDir() {
			continue
		}

		ext := strings.ToUpper(filepath.Ext(file.Name()))
		if ext == ".fa" || ext == ".fasta" {
			return file.Name(), nil
		}
	}

	return "", fmt.Errorf("failed: no input argument set and no fasta file found in %s", dir)
}

// parseFeatureInput turns the arguments to the features command into a
// CSV list of features in a single string.
func (p *inputParser) parseFeatureInput(args []string) (out string) {
	commaSeparated := false
	for _, a := range args {
		if strings.Contains(a, ",") {
			commaSeparated = true
		}
	}

	// if it's a features command, concatenate the arguments in case they're feature names
	// with 'repp features' the arguments are the feature names to use
	if commaSeparated {
		spacedSeparated := strings.Join(args, " ")
		splitByComma := strings.Split(spacedSeparated, ",")
		trimmedSeqs := []string{}
		for _, entry := range splitByComma {
			trimmedSeqs = append(trimmedSeqs, strings.TrimSpace(entry))
		}
		return strings.Join(trimmedSeqs, ",")
	}

	return strings.Join(args, " ")
}

// guessOutput gets an outpath path from an input path (if no output path is
// specified). It uses the same name as the input path to create an output.
func (p *inputParser) guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".output.json"
}

// parseDBs returns a list of absolute paths to BLAST databases.
func (p *inputParser) parseDBs(m *manifest, dbInput string) (dbs []DB, err error) {
	dbNames := p.parseCommaList(dbInput)

	if m.empty() {
		return nil, errors.New("no databases loaded. See 'repp add database'")
	}

	// if none filtered for, return all databases
	if len(dbNames) == 0 {
		for _, db := range m.DBs {
			dbs = append(dbs, db)
		}
		return
	}

	// filter for matching databases, error if none present
	for _, dbName := range dbNames {
		db, ok := m.DBs[dbName]
		if ok {
			dbs = append(dbs, db)
		} else {
			return nil, fmt.Errorf(
				"failed to find a DB with name: %s\n\tknown DBs: %s",
				dbName,
				strings.Join(m.GetNames(), ","),
			)
		}
	}

	return dbs, nil
}

// parseCommaList converts a comma separated list of strings into a list
// of strings for use elsewhere
func (p *inputParser) parseCommaList(commaList string) (newList []string) {
	comma := regexp.MustCompile(",")
	for _, item := range comma.Split(commaList, -1) {
		item = strings.Trim(item, " ,")
		if item == "" {
			continue
		}
		newList = append(newList, item)
	}

	return newList
}

// parseBackbone takes a backbone, referenced by its id, and an enzyme to cleave the
// backbone, and returns the linearized backbone as a Frag.
func (p *inputParser) parseBackbone(
	bbName string,
	enzymeNames []string,
	dbs []DB,
	c *config.Config,
) (f *Frag, backbone *Backbone, err error) {
	// if no backbone was specified, return an empty Frag
	if bbName == "" {
		return &Frag{}, &Backbone{}, nil
	}

	// confirm that the backbone exists in one of the dbs (or local fs) gather it as a Frag if it does
	bbFrag, err := queryDatabases(bbName, dbs)
	if err != nil {
		return &Frag{}, &Backbone{}, err
	}

	// try to digest the backbone with the enzyme
	if len(enzymeNames) == 0 {
		return &Frag{},
			&Backbone{},
			fmt.Errorf("backbone passed, %s, without an enzyme to digest it", bbName)
	}

	// gather the enzyme by name, err if it's unknown
	enzymes, err := p.getEnzymes(enzymeNames)
	if err != nil {
		return &Frag{}, &Backbone{}, err
	}

	if f, backbone, err = digest(bbFrag, enzymes); err != nil {
		return &Frag{}, &Backbone{}, err
	}

	return
}

// getEnzymes return the enzyme with the name passed. errors out if there is none.
func (p *inputParser) getEnzymes(enzymeNames []string) (enzymes []enzyme, err error) {
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

// getFilters takes an input string and returns a list of strings to run against matches
// when filtering out possible building fragments.
func (p *inputParser) getFilters(filterFlag string) []string {
	splitFunc := func(c rune) bool {
		return c == ' ' || c == ',' // space or comma separated
	}

	return strings.FieldsFunc(strings.ToUpper(filterFlag), splitFunc)
}

// read a FASTA file (by its path on local FS) to a slice of Fragments.
func read(path string, feature bool) (fragments []*Frag, err error) {
	if !filepath.IsAbs(path) {
		path, err = filepath.Abs(path)
		if err != nil {
			return nil, fmt.Errorf("failed to create path to input file: %s", err)
		}
	}

	dat, err := ioutil.ReadFile(path)
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
	id := idRegex.FindString(genbankSplit[0])

	if id == "" {
		return nil, fmt.Errorf("failed to parse locus from %s", path)
	}

	return []*Frag{
		{
			ID:  id,
			Seq: cleanedSeq,
		},
	}, nil
}
