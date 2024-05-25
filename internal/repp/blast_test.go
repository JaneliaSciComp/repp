package repp

import (
	"reflect"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

// test blast for circular genome with a left margin specified
func Test_BLAST_CircularGenomeWithLeftMargin(t *testing.T) {
	// create mock test fragment
	id := "test_target"
	seq := "GGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATGACCACCTTGATCTTCATGGTCTGGGTGCCCTCGTAGGGCTTGCCTTCGCCCTCGGATGTGCACTTGAAGTGGTGGTTGTTCACGGTGCCCTCCATGTACAGCTTCATGTGCATGTTCTCCTTGATCAGCTCGCTCATAGGTCCAGGGTTCTCCTCCACGTCTCCAGCCTGCTTCAGCAGGCTGAAGTTAGTAGCTCCGCTTCCGGATCCCCCGGGGAGCATGTCAAGGTCAAAATCGTCAAGAGCGTCAGCAGGCAGCATATCAAGGTCAAAGTCGTCAAGGGCATCGGCTGGGAgCATGTCTAAgTCAAAATCGTCAAGGGCGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGCTTACAATTTAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAA"
	leftMargin := 500

	// run blast
	matches, err := blast(id, seq, true, leftMargin, []DB{testDB}, []string{}, 10, false) // any match over 10 bp

	// check if it fails
	if err != nil {
		t.Errorf("failed to run BLAST: %v", err)
		return
	}

	// make sure matches are found
	if len(matches) < 1 {
		t.Error("failed to find any matches")
		return
	}

	for _, m := range matches {
		if m.queryStart < leftMargin {
			t.Errorf("match %v not expected in fragment matches", m)
		}
	}
}

// test the ability to find test fragments in a mock database
func Test_BLAST(t *testing.T) {
	// create mock test fragment
	id := "test_target"
	seq := "GGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATGACCACCTTGATCTTCATGGTCTGGGTGCCCTCGTAGGGCTTGCCTTCGCCCTCGGATGTGCACTTGAAGTGGTGGTTGTTCACGGTGCCCTCCATGTACAGCTTCATGTGCATGTTCTCCTTGATCAGCTCGCTCATAGGTCCAGGGTTCTCCTCCACGTCTCCAGCCTGCTTCAGCAGGCTGAAGTTAGTAGCTCCGCTTCCGGATCCCCCGGGGAGCATGTCAAGGTCAAAATCGTCAAGAGCGTCAGCAGGCAGCATATCAAGGTCAAAGTCGTCAAGGGCATCGGCTGGGAgCATGTCTAAgTCAAAATCGTCAAGGGCGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGCTTACAATTTAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAA"

	// run blast
	matches, err := blast(id, seq, true, 0, []DB{testDB}, []string{}, 10, false) // any match over 10 bp

	// check if it fails
	if err != nil {
		t.Errorf("failed to run BLAST: %v", err)
		return
	}

	// make sure matches are found
	if len(matches) < 1 {
		t.Error("failed to find any matches")
		return
	}

	matchesContain := func(targ match) {
		for _, m := range matches {
			if targ.entry == m.entry && targ.queryStart == m.queryStart && targ.queryEnd == m.queryEnd {
				return
			}
		}

		t.Errorf("failed to find match %v in fragment matches", targ)
	}

	matchesContain(match{
		entry:      "gnl|addgene|107006",
		queryStart: 0,
		queryEnd:   72,
	})
}

// test that we can filter out overlapping regions from blast results
// and those that are up against the edge of the fragment
func Test_cull(t *testing.T) {
	// test fragment with 3 matches that should be removed
	matches := []match{
		// shouldn't be removed
		{
			entry:      "m1",
			queryStart: 15,
			queryEnd:   19,
		},
		// should be removed because it fits within m3
		{
			entry:      "m2",
			queryStart: 29,
			queryEnd:   34,
		},
		// shouldn't be removed
		{
			entry:      "m3",
			queryStart: 29,
			queryEnd:   35,
		},
		// shouldn't be removed
		{
			entry:      "m4",
			queryStart: 31,
			queryEnd:   72,
		},
	}

	culledMatches := cull(matches, 3, 1)

	// make sure m2 has been removed
	for _, m := range culledMatches {
		if m.entry == "m2" {
			t.Error("m2 found in resulting matches, should have been removed")
		}
	}

	if len(culledMatches) != 3 {
		t.Errorf("%d filtered matches found on test fragment, 3 expected: %v", len(culledMatches), culledMatches)
	}

	if len(cull(matches, 3, 2)) != 4 {
		t.Errorf("%d filtered matches found on test fragment, 4 expected: %v", len(cull(matches, 3, 2)), cull(matches, 3, 2))
	}
}

func Test_isMismatch(t *testing.T) {
	c := config.New()
	c.PcrPrimerMaxOfftargetTm = 40.0

	type args struct {
		sequence string
		match    match
	}
	tests := []struct {
		name string
		args args
		want bool
	}{
		{
			"find mismatch",
			args{
				sequence: "gtccgcgtcgtcgtcat",
				match: match{
					seq:                 "atgacgacgacgcggac",
					queryRevCompMatch:   false,
					subjectRevCompMatch: true,
				},
			},
			true,
		},
		// this is only ~15deg
		{
			"no false positive mismatch",
			args{
				sequence: "gtccgcgtcgtcgtcat",
				match: match{
					seq:                 "acgacgacgac",
					queryRevCompMatch:   false,
					subjectRevCompMatch: true,
				},
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := isMismatch(tt.args.sequence, tt.args.match, c); got != tt.want {
				t.Errorf("isMismatch() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parentMismatch(t *testing.T) {
	conf := config.New()
	conf.PcrPrimerMaxOfftargetTm = 35.0

	type args struct {
		primer string
		parent string
	}
	tests := []struct {
		name         string
		args         args
		wantMismatch bool
		wantMatch    match
		wantErr      bool
	}{
		{
			"avoids false positive",
			args{
				"GTTGGAGTCCACGTTCTTT",
				"gnl|addgene|107006",
			},
			false,
			match{},
			false,
		},
		// I intentionally added another off-target seq to 107006, AGTATAGTAGGTAGTCATTCTT
		{
			"finds mismatch",
			args{
				"AGTATAGGATAGGTAGTCATTCTT",
				"gnl|addgene|107006",
			},
			true,
			match{
				entry:       "addgene:107006",
				uniqueID:    "addgene:107006-0",
				seq:         "AGTATAGTAGGTAGTCATTCTT",
				querySeq:    "AGTATAGGATAGGTAGTCATTCTT",
				queryStart:  0,
				queryEnd:    23,
				mismatching: 2,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			mismatchResult := parentMismatch([]Primer{{Seq: tt.args.primer}}, tt.args.parent, testDB, conf)
			gotMismatch := mismatchResult.wasMismatch
			gotMatch := mismatchResult.m
			err := mismatchResult.err

			if (err != nil) != tt.wantErr {
				t.Errorf("parentMismatch() error = %+v, wantErr %+v", err, tt.wantErr)
				return
			}
			if gotMismatch != tt.wantMismatch {
				t.Errorf("parentMismatch() gotMismatch = %+v, want %+v", gotMismatch, tt.wantMismatch)
			}

			// have to mutate the fields not included in expected set
			gotMatch.circular = false
			gotMatch.title = ""
			gotMatch.subjectStart = 0
			gotMatch.subjectEnd = 0
			gotMatch.queryRevCompMatch = false
			gotMatch.subjectRevCompMatch = false

			if !reflect.DeepEqual(gotMatch, tt.wantMatch) {
				t.Errorf("parentMismatch() gotMatch = %+v, want %+v", gotMatch, tt.wantMatch)
			}
		})
	}
}

func Test_queryDatabases(t *testing.T) {
	type args struct {
		entry string
		dbs   []DB
	}
	tests := []struct {
		name    string
		args    args
		wantF   Frag
		wantErr bool
	}{
		{
			name: "query 85039.2",
			args: args{
				entry: "addgene:85039.2",
				dbs:   []DB{testDB},
			},
			wantF: Frag{
				ID: "addgene:85039.2",
				db: testDB,
			},
			wantErr: false,
		},
		{
			name: "fail at nonsense frag",
			args: args{
				entry: "jahf9a8f9",
				dbs:   []DB{testDB},
			},
			wantF:   Frag{},
			wantErr: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotF, err := queryDatabases(tt.args.entry, tt.args.dbs)
			if (err != nil) != tt.wantErr {
				t.Errorf("queryDatabases() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if gotF.ID != tt.wantF.ID {
				t.Errorf("queryDatabases().ID = %v, want %v", gotF.ID, tt.wantF.ID)
			}
			if gotF.db != tt.wantF.db {
				t.Errorf("queryDatabases().DB = %v, want %v", gotF.db, tt.wantF.db)
			}
		})
	}
}

func Test_blastdbcmd(t *testing.T) {
	type args struct {
		entry string
		db    DB
	}
	tests := []struct {
		name    string
		args    args
		wantErr bool
	}{
		{
			name: "find 107006",
			args: args{
				entry: "gnl|addgene|107006",
				db:    testDB,
			},
			wantErr: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, _, err := blastdbcmd(tt.args.entry, tt.args.db)
			if (err != nil) != tt.wantErr {
				t.Errorf("blastdbcmd() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
		})
	}
}
