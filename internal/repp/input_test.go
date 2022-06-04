package repp

import (
	"path"
	"reflect"
	"testing"
)

func Test_inputParser_parseOut(t *testing.T) {
	parser := inputParser{}

	type args struct {
		in string
	}
	tests := []struct {
		name    string
		args    args
		wantOut string
	}{
		{
			"parse relative path to neighboring output path",
			args{
				in: "./test_file.fa",
			},
			"./test_file.output.json",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotOut := parser.guessOutput(tt.args.in); gotOut != tt.wantOut {
				t.Errorf("parseOut() = %v, want %v", gotOut, tt.wantOut)
			}
		})
	}
}

func Test_inputParser_getFilters(t *testing.T) {
	type args struct {
		filterFlag string
	}
	tests := []struct {
		name string
		p    *inputParser
		args args
		want []string
	}{
		{
			"biobrick separated from year by commas",
			&inputParser{},
			args{
				filterFlag: "tests,BBa_k222000,2004",
			},
			[]string{"TESTS", "BBA_K222000", "2004"},
		},
		{
			"single year",
			&inputParser{},
			args{
				filterFlag: "2004",
			},
			[]string{"2004"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &inputParser{}
			if got := p.getFilters(tt.args.filterFlag); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("inputParser.getFilters() = %v, want %v", got, tt.want)
			}
		})
	}
}

// Test reading of a FASTA file
func Test_read(t *testing.T) {
	type fileRead struct {
		name         string
		file         string
		fragCount    int
		readFeatures bool
	}

	files := []fileRead{
		{
			"113726(circular)",
			path.Join("..", "..", "test", "input", "113726(circular).parent"),
			1,
			false,
		},
		{
			"multi.fasta",
			path.Join("..", "..", "test", "input", "multi.fasta"),
			5,
			false,
		},
		{
			"genbank sequence",
			path.Join("..", "..", "test", "input", "genbank.gb"),
			1,
			false,
		},
		{
			"genbank features",
			path.Join("..", "..", "test", "input", "genbank.gb"),
			66,
			true,
		},
	}

	for _, f := range files {
		fragments, err := read(f.file, f.readFeatures)

		if err != nil {
			t.Error(err)
		}

		if len(fragments) != f.fragCount {
			t.Errorf("failed to load fragments, len=%d, expected=%d", len(fragments), f.fragCount)
		}

		for _, f := range fragments {
			// ensure we got an ID
			if len(f.ID) < 1 {
				t.Error("failed to load an ID for a Frag from FASTA")
			}

			// ensure we got a Seq
			if len(f.Seq) < 1 {
				t.Errorf("failed to parse a sequence for Frag %s", f.ID)
			}
		}
	}
}

func Test_inputParser_parseDBs(t *testing.T) {
	testDB := DB{
		Path: "/tmp/test.fa",
	}

	type args struct {
		manifest *manifest
		dbInput string
	}
	tests := []struct {
		name    string
		p       *inputParser
		args    args
		wantDbs []DB
		wantErr bool
	}{
		{
			name: "return all by default",
			p: &inputParser{},
			args: args{
				manifest: &manifest{
					DBs: map[string]DB{
						"test": testDB,
					},
				},
				dbInput: "",
			},
			wantDbs: []DB{testDB},
			wantErr: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &inputParser{}
			gotDbs, err := p.parseDBs(tt.args.manifest, tt.args.dbInput)
			if (err != nil) != tt.wantErr {
				t.Errorf("inputParser.parseDBs() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(gotDbs, tt.wantDbs) {
				t.Errorf("inputParser.parseDBs() = %v, want %v", gotDbs, tt.wantDbs)
			}
		})
	}
}
