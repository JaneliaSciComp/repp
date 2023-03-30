package repp

import (
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/Lattice-Automation/repp/internal/config"
	"go.uber.org/multierr"
)

// match is a blast "hit" in the blastdb.
type match struct {
	// entry of the matched building fragment in the database
	entry string

	// unique id for the match (entry name + start index % seqL). also used by fragments
	uniqueID string

	// seq of the match on the subject match
	querySeq string

	// queryStart of the queried seq match (0-indexed)
	queryStart int

	// queryEnd of the queried seq match (0-indexed)
	queryEnd int

	// sequence of the match in the subject sequence
	seq string

	// start index of the match on the subject fragment (entry) (0-indexed)
	subjectStart int

	// end index of the match on the subject fragment (entry) (0-indexed)
	subjectEnd int

	// the db that was BLASTed against (used later for checking off-targets in parents)
	db DB

	// titles from the db. eg: year it was created
	title string

	// circular if it's a circular fragment in the db (plasmid, plasmid, etc)
	circular bool

	// mismatching number of bps in the match (for primer off-targets)
	mismatching int

	// forward if the match is along the sequence strand versus the reverse complement strand
	forward bool
}

// length returns the length of the match on the queried fragment.
func (m *match) length() int {
	queryLength := m.queryEnd - m.queryStart + 1
	subjectLength := m.subjectEnd - m.subjectStart + 1

	if queryLength > subjectLength {
		return queryLength
	}

	return subjectLength
}

// mismatchResults are the results of a seqMismatch check. saved between runs for perf
type mismatchResult struct {
	wasMismatch bool
	m           match
	err         error
}

// blastExec is a small utility object for executing BLAST.
type blastExec struct {
	// the name of the query
	name string

	// the sequence of the query
	seq string

	// whether to circularize the queries sequence in the input file
	circular bool

	// the database we're BLASTing against
	db DB

	// the input BLAST file
	in *os.File

	// the output BLAST file
	out *os.File

	// optional path to a FASTA file with a subject FASTA sequence
	subject string

	// the percentage identity for BLAST queries
	identity int

	// the expect value of a BLAST query (defaults to 10)
	evalue int
}

// input creates an input query file (FASTA) for blastn.
func (b *blastExec) input() error {
	// create the query sequence file.
	var querySeq string

	if b.circular {
		// if circular, add the sequence to itself because it's circular
		// and we want to find matches across the zero-index
		querySeq = b.seq + b.seq
	} else {
		querySeq = b.seq
	}

	_, err := b.in.WriteString(fmt.Sprintf(">%s\n%s\n", b.name, querySeq))

	return err
}

// run calls the external blastn binary on the input file.
func (b *blastExec) run() (err error) {
	threads := runtime.NumCPU() - 1
	if threads < 1 {
		threads = 1
	}

	rlog.Infof("Query %s against %s -> %s\n", b.in.Name(),
		b.db.Path, b.out.Name())

	flags := []string{
		"-task", "blastn",
		"-db", b.db.Path,
		"-query", b.in.Name(),
		"-out", b.out.Name(),
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch gaps stitle",
		"-perc_identity", fmt.Sprintf("%d", b.identity),
		// "-culling_limit", "50",
		// "-max_target_seqs", "5000",
		"-num_threads", strconv.Itoa(threads),
	}

	if b.identity > 99 {
		flags = append(flags,
			"-reward", "1", // most negative penalty I could find in blast_stat.c
			"-penalty", "-5", // needed because mismatches were being included in the end of pSB1A3 matches
			"-gapopen", "6",
			"-gapextend", "6",
		)
	} else if b.identity >= 98 {
		// Table D1 https://www.ncbi.nlm.nih.gov/books/NBK279684/
		flags = append(flags,
			"-reward", "1",
			"-penalty", "-3",
			"-gapopen", "3",
			"-gapextend", "3",
		)
	} else if b.identity >= 90 {
		// Table D1 https://www.ncbi.nlm.nih.gov/books/NBK279684/
		flags = append(flags,
			"-reward", "1",
			"-penalty", "-2",
			"-gapopen", "1",
			"-gapextend", "2",
		)
	} else {
		// Table D1 https://www.ncbi.nlm.nih.gov/books/NBK279684/
		flags = append(flags,
			"-reward", "1",
			"-penalty", "-1",
			"-gapopen", "1",
			"-gapextend", "2",
		)
	}

	if b.evalue != 0 {
		flags = append(flags, "-evalue", strconv.Itoa(b.evalue))
	} else if b.identity < 90 {
		flags = append(flags, "-evalue", "5000")
	} else if b.identity < 98 {
		flags = append(flags, "-evalue", "1000")
	}

	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command("blastn", flags...)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to execute blastn against %s: %v: %s", b.db.Name, err, string(output))
	}

	return
}

// parse reads the output of blastn into matches.
func (b *blastExec) parse(filters []string) (matches []match, err error) {
	// read in the results
	file, err := os.ReadFile(b.out.Name())
	if err != nil {
		return
	}
	fileS := string(file)

	fullQuery := b.seq + b.seq
	identityThreshold := float64(b.identity)/100.0 - 0.0001

	// read it into Matches
	var ms []match
	for _, line := range strings.Split(fileS, "\n") {
		// comment lines start with a #
		if strings.HasPrefix(line, "#") {
			continue
		}

		// split on white space
		cols := strings.Fields(line)
		if len(cols) < 6 {
			continue
		}

		entry := strings.Replace(cols[0], ">", "", -1)
		queryStart, _ := strconv.Atoi(cols[1])
		queryEnd, _ := strconv.Atoi(cols[2])
		subjectStart, _ := strconv.Atoi(cols[3])
		subjectEnd, _ := strconv.Atoi(cols[4])
		seq := cols[5]                          // subject sequence
		mismatching, _ := strconv.Atoi(cols[6]) // mismatch count
		gaps, _ := strconv.Atoi(cols[7])        // gap count
		titles := cols[8]                       // salltitles, eg: "fwd-terminator-2011"
		forward := true

		// check whether the mismatch ratio is less than the set limit
		matchRatio := float64(len(seq)-(mismatching+gaps)) / float64(len(seq))
		if matchRatio < identityThreshold {
			continue
		}

		seq = strings.Replace(seq, "-", "", -1) // remove gap markers
		queryStart--                            // convert from 1-based to 0-based
		queryEnd--
		subjectStart--
		subjectEnd--

		// bug where titles are being included in the entry
		entryCols := strings.Fields(entry)
		if len(entryCols) > 1 {
			entry = entryCols[0]
			titles = entryCols[1] + titles
		}

		// flip if blast is reading right to left
		if queryStart > queryEnd {
			queryStart, queryEnd = queryEnd, queryStart
			forward = !forward
		}
		if subjectStart > subjectEnd {
			subjectStart, subjectEnd = subjectEnd, subjectStart
			forward = !forward
		}

		// filter on titles
		matchesFilter := false
		titles += entry
		titles = strings.ToUpper(titles)
		for _, f := range filters {
			if strings.Contains(titles, f) {
				matchesFilter = true
				break
			}
		}
		if matchesFilter {
			continue // has been filtered out because of the "exclude" CLI flag
		}

		// get a unique identifier to distinguish this match/fragment from the others
		uniqueID := entry + strconv.Itoa(queryStart%len(b.seq))

		// gather the query sequence
		querySeq := fullQuery[queryStart : queryEnd+1]

		// create and append the new match
		ms = append(ms, match{
			entry:        entry,
			uniqueID:     uniqueID,
			querySeq:     querySeq,
			queryStart:   queryStart,
			queryEnd:     queryEnd,
			seq:          seq,
			subjectStart: subjectStart,
			subjectEnd:   subjectEnd,
			circular:     strings.Contains(entry+titles, "CIRCULAR"),
			mismatching:  mismatching + gaps,
			db:           b.db,
			title:        titles,
			forward:      forward,
		})
	}

	return ms, nil
}

// runs blast on the query file against another subject file (rather than blastdb)
func (b *blastExec) runAgainst() (err error) {
	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		"blastn",
		"-task", "blastn",
		"-query", b.in.Name(),
		"-subject", b.subject,
		"-out", b.out.Name(),
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch gaps stitle",
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to execute blastn against %s: %v: %s", b.subject, err, string(output))
	}
	return
}

func (b *blastExec) close() error {
	var result error
	result = multierr.Append(result, os.Remove(b.in.Name()))
	result = multierr.Append(result, os.Remove(b.out.Name()))
	return result

}

// blast the seq against all dbs and acculate matches.
func blast(
	name, seq string,
	circular bool,
	dbs []DB,
	filters []string,
	identity int,
) ([]match, error) {
	matches := []match{}
	for _, db := range dbs {
		in, err := os.CreateTemp("", "blast-in-*")
		if err != nil {
			return nil, err
		}

		out, err := os.CreateTemp("", "blast-out-*")
		if err != nil {
			return nil, err
		}

		b := &blastExec{
			name:     name,
			seq:      seq,
			circular: circular,
			db:       db,
			in:       in,
			out:      out,
			identity: identity,
		}
		defer b.close()

		// make sure the db exists
		if _, err := os.Stat(db.Path); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find a BLAST database at %s", db.Path)
		}

		// create the input file
		if err := b.input(); err != nil {
			return nil, fmt.Errorf("failed to write a BLAST input file at %s: %v", b.in.Name(), err)
		}

		// execute BLAST
		if err := b.run(); err != nil {
			return nil, fmt.Errorf("failed executing BLAST: %v", err)
		}

		// parse the output file to Matches against the Frag
		rlog.Infof("Parse filters %+q", filters)
		dbMatches, err := b.parse(filters)
		if err != nil {
			return nil, fmt.Errorf("failed to parse BLAST output: %v", err)
		}

		// add these matches against the growing list of matches
		matches = append(matches, dbMatches...)
	}

	return matches, nil
}

// blastAgainst runs against a pre-made subject database
func blastAgainst(
	name, seq, subject string,
	circular bool,
	identity int,
) (matches []match, err error) {
	in, err := os.CreateTemp("", "blast-in-*")
	if err != nil {
		return nil, err
	}

	out, err := os.CreateTemp("", "blast-out-*")
	if err != nil {
		return nil, err
	}

	b := &blastExec{
		name:     name,
		seq:      seq,
		circular: circular,
		subject:  subject,
		in:       in,
		out:      out,
		identity: identity,
	}
	defer b.close()

	// make sure the subject file exists
	if _, err := os.Stat(subject); os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to find a BLAST subject at %s", subject)
	}

	// create the input file
	if err := b.input(); err != nil {
		return nil, fmt.Errorf("failed to write a BLAST input file at %s: %v", b.in.Name(), err)
	}

	// execute BLAST
	if err := b.runAgainst(); err != nil {
		return nil, fmt.Errorf("failed executing BLAST: %v", err)
	}

	// parse the output file to Matches against the Frag
	if matches, err = b.parse([]string{}); err != nil {
		return nil, fmt.Errorf("failed to parse BLAST output: %v", err)
	}

	return matches, nil
}

// cull removes matches that are engulfed in others
//
// culling fragment matches means removing those that are completely
// self-contained in other fragments
// if limit == 1 the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
func cull(matches []match, targetLength, minSize, limit int) (culled []match) {
	// remove fragments that are shorter the minimum cut off size
	// propertize by source because this isn't smart enough to propertize based on
	// cost as well. Ie we want to avoid propertizing a small fragment enclosed in a larger
	// but far more expensive fragment (and instead calculate that later)
	groupedMatches := map[string][]match{}
	for _, m := range matches {
		if minSize > 0 && m.length() < minSize {
			continue // too short
		}
		groupedMatches[m.db.Path] = append(groupedMatches[m.db.Path], m)
	}

	// create culled matches (non-self contained)
	for _, group := range groupedMatches {
		culled = append(culled, properize(group, limit)...)
	}

	// because we culled the matches, we may have removed a match from the
	// start or the end. right now, a match showing up twice in the plasmid
	// is how we circularize, so have to add back matches to the start or end
	matchCount := make(map[string]int)
	for _, m := range culled {
		matchCount[m.uniqueID]++
	}

	// sort again now that we added copied matches
	sortMatches(culled)
	return culled
}

// properize remove matches that are entirely contained within others
func properize(matches []match, limit int) []match {
	sortMatches(matches)

	// only include those that aren't encompassed by the one before it
	culled := []match{}
	for _, m := range matches {
		check := len(culled) - limit
		if check < 0 || m.queryEnd > culled[check].queryEnd {
			culled = append(culled, m)
		}
	}

	return culled
}

// sortMatches sorts matches by their start index
// for fragments with equivelant starting indexes, put the larger one first
func sortMatches(matches []match) {
	sort.Slice(matches, func(i, j int) bool {
		if matches[i].queryStart != matches[j].queryStart {
			return matches[i].queryStart < matches[j].queryStart
		} else if matches[i].length() != matches[j].length() {
			return matches[i].length() > matches[j].length()
		} else if matches[i].circular && !matches[j].circular {
			return true
		} else if !matches[i].circular && matches[j].circular {
			return false
		}
		return matches[i].entry > matches[j].entry
	})
}

// queryDatabases is for finding a fragment/plasmid with the entry name in one of the dbs
func queryDatabases(entry string, dbs []DB) (f *Frag, err error) {
	// first try to get the entry out of a local file
	if frags, err := read(entry, false); err == nil && len(frags) > 0 {
		return frags[0], nil // it was a local file
	}

	// channel that returns filename to an output result from blastdbcmd
	outFileCh := make(chan string, len(dbs))
	dbSourceCh := make(chan DB, len(dbs))

	// move through each db and see if it contains the entry
	for _, db := range dbs {
		go func(db DB) {
			// if outFile is defined here we managed to query the entry from the db
			outFile, _, err := blastdbcmd(entry, db)
			if err == nil && outFile != nil {
				outFileCh <- outFile.Name() // "" if not found
				dbSourceCh <- db
			} else {
				outFileCh <- ""
				dbSourceCh <- db
			}
		}(db)
	}

	// try and return a read fragment file
	for i := 0; i < len(dbs); i++ {
		outFile := <-outFileCh
		dbSource := <-dbSourceCh
		if outFile == "" {
			continue // failed to query from this DB
		}
		defer os.Remove(outFile)

		if frags, err := read(outFile, false); err == nil {
			targetFrag := frags[0]

			// fix the ID, don't want titles in the ID (bug)
			idSplit := strings.Fields(targetFrag.ID)
			if len(idSplit) > 1 {
				targetFrag.ID = idSplit[0]
			}

			targetFrag.db = dbSource
			return targetFrag, nil
		}

		return &Frag{}, err
	}

	close(outFileCh)
	close(dbSourceCh)

	return &Frag{}, fmt.Errorf("failed to find frag %s in any of: %s", entry, strings.Join(dbNames(dbs), ","))
}

// seqMismatch queries for any mismatching primer locations in the parent sequence
// unlike parentMismatch, it doesn't first find the parent fragment from the db it came from
// the sequence is passed directly as parentSeq
func seqMismatch(primers []Primer, parentID, parentSeq string, conf *config.Config) mismatchResult {
	parentFile, err := os.CreateTemp("", "parent-*")
	if err != nil {
		return mismatchResult{false, match{}, err}
	}
	defer os.Remove(parentFile.Name())

	if parentID == "" {
		parentID = "parent"
	}
	inContent := fmt.Sprintf(">%s\n%s\n", parentID, parentSeq)
	if _, err = parentFile.WriteString(inContent); err != nil {
		return mismatchResult{false, match{}, fmt.Errorf("failed to write primer sequence to query FASTA file: %v", err)}
	}

	// check each primer for mismatches
	for _, primer := range primers {
		wasMismatch, m, err := mismatch(primer.Seq, parentFile, conf)
		if wasMismatch || err != nil {
			return mismatchResult{wasMismatch, m, err}
		}
	}

	return mismatchResult{false, match{}, nil}
}

// parentMismatch both searches for a the parent fragment in its source DB and queries for
// any mismatches in the seq before returning
func parentMismatch(primers []Primer, parent string, db DB, conf *config.Config) mismatchResult {
	// try and query for the parent in the source DB and write to a file
	parentFile, parentSeq, err := blastdbcmd(parent, db)

	// ugly check here for whether we just failed to get the parent entry from a db
	// which isn't a huge deal (shouldn't be flagged as a mismatch)
	// this is similar to what io.IsNotExist does
	if err != nil {
		if strings.Contains(err.Error(), "failed to query") {
			stderr.Println(err) // just write the error
			// TODO: if we fail to find the parent, query the fullSeq as it was sent
			return mismatchResult{false, match{}, nil}
		}
		return mismatchResult{false, match{}, err}
	}

	// check each primer for mismatches
	if parentFile.Name() != "" {
		defer os.Remove(parentFile.Name())

		for i, primer := range primers {
			// confirm that the 3' end of the primer is in the parent seq
			primerEnd := primer.Seq[len(primer.Seq)-10:]
			if !strings.Contains(parentSeq, primerEnd) && !strings.Contains(parentSeq, reverseComplement(primerEnd)) {
				dir := "FWD"
				if i > 0 {
					dir = "REV"
				}
				return mismatchResult{false, match{}, fmt.Errorf("does not contain end of %s primer: %s", dir, primerEnd)}
			}

			// check for a mismatch in the parent sequence
			wasMismatch, m, err := mismatch(primer.Seq, parentFile, conf)
			if wasMismatch || err != nil {
				return mismatchResult{wasMismatch, m, err}
			}
		}
	}

	return mismatchResult{false, match{}, err}
}

// blastdbcmd queries a fragment/plasmid by its FASTA entry name (entry) and writes the
// results to a temporary file (to be BLAST'ed against)
//
// entry here is the ID that's associated with the fragment in its source DB (db)
func blastdbcmd(entry string, db DB) (output *os.File, parentSeq string, err error) {
	// path to the entry batch file to hold the entry accession
	entryFile, err := os.CreateTemp("", "blastcmd-in-*")
	if err != nil {
		return nil, "", err
	}
	defer os.Remove(entryFile.Name())

	// path to the output sequence file from querying the entry's sequence from the BLAST db
	output, err = os.CreateTemp("", "blastcmd-out-*")
	if err != nil {
		return nil, "", err
	}

	// write entry to file
	// this was a 2-day issue I couldn't resolve...
	// I was using the "-entry" flag on exec.Command, but have since
	// switched to the simpler -entry_batch command (on a file) that resolves the issue
	if _, err := entryFile.WriteString(entry); err != nil {
		return nil, "", fmt.Errorf("failed to write blastdbcmd entry file at %s: %v", entryFile.Name(), err)
	}

	// make a blastdbcmd command (for querying a DB, very different from blastn)
	queryCmd := exec.Command(
		"blastdbcmd",
		"-db", db.Path,
		"-dbtype", "nucl",
		"-entry_batch", entryFile.Name(),
		"-out", output.Name(),
		"-outfmt", "%f ", // fasta format
	)

	// execute
	if _, err := queryCmd.CombinedOutput(); err != nil {
		return nil, "", fmt.Errorf("warning: failed to query %s from %s db\n\t%s", entry, db.Name, err.Error())
	}

	// read in the results as fragments. set their sequence to the full one returned from blastdbcmd
	fragments, err := read(output.Name(), false)
	if err == nil && len(fragments) >= 1 {
		for _, f := range fragments {
			f.fullSeq = f.Seq // set fullSeq, faster to check for primer off-targets later
			return output, f.Seq, nil
		}
	}

	return nil, "", fmt.Errorf("warning: failed to query %s from %s db", entry, db.Name)
}

// mismatch finds mismatching sequences between the query sequence and
// the parent sequence (in the parent file)
//
// The fragment to query against is stored in parentFile
func mismatch(primer string, parentFile *os.File, c *config.Config) (wasMismatch bool, m match, err error) {
	// path to the entry batch file to hold the entry accession
	in, err := os.CreateTemp("", "primer3-in-*")
	if err != nil {
		return false, match{}, err
	}

	// path to the output sequence file from querying the entry's sequence from the BLAST db
	out, err := os.CreateTemp("", "primer3-out-*")
	if err != nil {
		return false, match{}, err
	}

	// create input file
	inContent := fmt.Sprintf(">primer\n%s\n", primer)
	if _, err = in.WriteString(inContent); err != nil {
		return false, m, fmt.Errorf("failed to write primer sequence to query FASTA file: %v", err)
	}

	// BLAST the query sequence against the parentFile sequence
	b := &blastExec{
		in:       in,
		out:      out,
		subject:  parentFile.Name(),
		seq:      primer,
		identity: 65,    // see Primer-BLAST https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3412702/
		evalue:   30000, // see Primer-BLAST
	}
	defer b.close()

	// execute BLAST
	if err = b.runAgainst(); err != nil {
		return false, m, fmt.Errorf("failed to run blast against parent: %v", err)
	}

	// get the BLAST matches
	matches, err := b.parse([]string{})
	if err != nil {
		return false, match{}, fmt.Errorf("failed to parse matches from %s: %v", out.Name(), err)
	}

	// parse the results and check whether any are cause for concern (by Tm)
	primerCount := 1 // number of times we expect to see the primer itself
	parentFileContents, err := os.ReadFile(parentFile.Name())
	if err != nil {
		return false, match{}, err
	}

	if strings.Contains(string(parentFileContents), "circular") {
		// if the match is against a circular fragment, we expect to see the primer's binding location
		// twice because circular fragments' sequences are doubled in the DBs
		primerCount++
	}

	for _, m := range matches {
		if isMismatch(primer, m, c) {
			primerCount--
		}

		if primerCount < 0 {
			return true, m, nil
		}
	}

	return false, match{}, nil
}

// isMismatch returns whether the match constitutes a mismatch
// between it and the would be primer sequence
//
// estimate the ntthal and check against the max offtarget tm
// from the settings
func isMismatch(primer string, m match, c *config.Config) bool {
	// we want the reverse complement of one to the other
	ectopic := m.seq
	if m.forward {
		ectopic = reverseComplement(ectopic)
	}

	ntthalCmd := exec.Command(
		"ntthal",
		"-a", "END1", // end of primer sequence
		"-s1", primer,
		"-s2", ectopic,
		"-path", config.Primer3Config,
		"-r", // temperature only
	)

	ntthalOut, err := ntthalCmd.CombinedOutput()
	if err != nil {
		stderr.Printf("failed to execute ntthal: %s", strings.Join(ntthalCmd.Args, ","))
		return true
	}

	ntthalOutString := string(ntthalOut)
	temp, err := strconv.ParseFloat(strings.TrimSpace(ntthalOutString), 64)
	if err != nil {
		rlog.Fatal(err)
	}

	return temp > c.PcrPrimerMaxOfftargetTm
}

// makeblastdb runs makeblastdb against a FASTA file.
func makeblastdb(fullDbPath string) error {
	rlog.Infof("Make BlastDB %s\n", fullDbPath)
	cleanblastdb(fullDbPath)

	cmd := exec.Command("makeblastdb",
		"-dbtype", "nucl",
		"-in", fullDbPath,
		"-parse_seqids",
	)

	rlog.Debugf("Run: %v", cmd.Args)
	if stdout, err := cmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to makeblastdb: %s %w", string(stdout), err)
	}
	return nil
}

func cleanblastdb(fullDbPath string) {
	blastdbExts := []string{
		".nhr", ".nos", ".nto", ".nin", ".not",
		".njs", ".nsq", ".ndb", ".nog", ".ntf",
	}

	for _, ext := range blastdbExts {
		dbFilename := fullDbPath + ext
		rlog.Debugf("Delete %s", dbFilename)
		if err := os.Remove(dbFilename); err != nil {
			rlog.Debugf("Error: %v", err)
		}
	}
}
