package repp

import (
	"reflect"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

// func Test_sequence_e2e(test *testing.T) {
// 	c := config.New()

// 	type testFlags struct {
// 		in       string
// 		out      string
// 		backbone string
// 		enzymes  []string
// 		filters  string
// 		dbs      []DB
// 	}

// 	tests := []testFlags{
// 		{
// 			path.Join("..", "..", "test", "input", "backbone.fa"),
// 			path.Join("..", "..", "test", "output", "backbone.json"),
// 			"pSB1A3",
// 			[]string{"PstI"},
// 			"2018,2019",
// 			[]DB{testDB},
// 		},
// 		{
// 			path.Join("..", "..", "test", "input", "BBa_K2224001.fa"),
// 			path.Join("..", "..", "test", "output", "BBa_K2224001.json"),
// 			"pSB1A3",
// 			[]string{"PstI"},
// 			"2017,2018,2019",
// 			[]DB{testDB},
// 		},
// 	}

// 	for _, t := range tests {
// 		sols := Sequence(NewFlags(t.in, t.out, t.backbone, t.filters, t.enzymes, t.dbs))

// 		if len(sols) < 1 {
// 			test.Errorf("no solutions for %s", t.in)
// 		}

// 		for _, s := range sols {
// 			e := validateJunctions(s, c)
// 			if e != nil {
// 				test.Logf("failed making %s\n", t.in)
// 				test.Error(e)
// 			}
// 		}
// 	}
// }

// func Test_features(t *testing.T) {
// 	test1, conf := NewFlags(
// 		"p10 promoter, mEGFP, T7 terminator",
// 		filepath.Join("..", "..", "test", "output", "features.json"),
// 		"pSB1A3",
// 		"",
// 		[]string{"EcoRI"},
// 		[]DB{},
// 	)

// 	test2, _ := NewFlags(
// 		"BBa_R0062,BBa_B0034,BBa_C0040,BBa_B0010,BBa_B0012",
// 		filepath.Join("..", "..", "test", "output", "igem.features.json"),
// 		"pSB1C3",
// 		"",
// 		[]string{"PstI", "EcoRI"},
// 		[]DB{},
// 	)

// 	type args struct {
// 		flags *Flags
// 		conf  *config.Config
// 	}
// 	tests := []struct {
// 		name string
// 		args args
// 	}{
// 		{
// 			"test end to end features creation",
// 			args{
// 				flags: test1,
// 				conf:  conf,
// 			},
// 		},
// 		{
// 			"test end to end features creation using iGEM parts",
// 			args{
// 				flags: test2,
// 				conf:  conf,
// 			},
// 		},
// 	}
// 	for _, tt := range tests {
// 		t.Run(tt.name, func(t *testing.T) {
// 			sols := Features(tt.args.flags, tt.args.conf)

// 			if len(sols) < 1 {
// 				t.Failed()
// 			}

// 			for _, s := range sols {
// 				e := validateJunctions(s, conf)
// 				if e != nil {
// 					t.Error(e)
// 				}
// 			}
// 		})
// 	}
// }

func Test_fragments(t *testing.T) {
	c := config.New()
	c.PcrMinLength = 10
	c.FragmentsMinHomology = 8
	c.FragmentsMaxHomology = 20

	type args struct {
		inputFragments []*Frag
		conf           *config.Config
	}
	tests := []struct {
		name              string
		args              args
		wantTargetPlasmid *Frag
		wantFragments     []*Frag
	}{
		{
			"fragments with linear overlap",
			args{
				[]*Frag{
					{
						Seq:  "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
						conf: c,
					},
					{
						Seq:  "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
						conf: c,
					},
					{
						Seq:  "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
						conf: c,
					},
				},
				c,
			},
			&Frag{
				Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCGACTGATCACTAGCATCGACTAGCTAGAACTGATCTAG",
			},
			[]*Frag{
				{
					Seq:      "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					fragType: linear,
				},
				{
					Seq:      "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					fragType: linear,
				},
				{
					Seq:      "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					fragType: linear,
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotTargetPlasmid, gotFragments := fragments(tt.args.inputFragments, tt.args.conf)

			if !reflect.DeepEqual(gotTargetPlasmid.Seq, tt.wantTargetPlasmid.Seq) {
				t.Errorf("fragments() gotTargetPlasmid = %v, want %v", gotTargetPlasmid, tt.wantTargetPlasmid)
			}

			if len(gotFragments) != len(tt.wantFragments) {
				t.Errorf("fragments() got %d fragments, expected %d", len(gotFragments), len(tt.wantFragments))
				return
			}

			for i, wantF := range tt.wantFragments {
				if wantF.Seq != gotFragments[i].Seq {
					t.Errorf("fragments() gotFragment.Seq = %v, want %v", gotFragments[i].Seq, wantF.Seq)
				}

				if wantF.fragType != gotFragments[i].fragType {
					t.Errorf("fragments() gotFragment.Type = %v, want %v", gotFragments[i].fragType, wantF.fragType)
				}
			}
		})
	}
}

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
// func Test_plasmid_single_plasmid(t *testing.T) {
// 	dir := t.TempDir()
// 	output, err := ioutil.TempFile(dir, "")
// 	if err != nil {
// 		t.Error(err)
// 	}

// 	fs, c := NewFlags(
// 		path.Join("..", "..", "test", "input", "109049.addgene.fa"),
// 		output.Name(),
// 		"",
// 		"",
// 		[]string{},
// 		[]DB{testDB},
// 	)

// 	assemblies := Sequence(fs, c)

// 	if !strings.Contains(assemblies[0][0].ID, "109049") {
// 		t.Fatal("failed to use 109049 to build the plasmid")
// 	}
// }
