package repp

import (
	"path"
	"reflect"
	"strings"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

type mockAssemblyParams struct {
	assemblyParamsImpl
}

func (ap mockAssemblyParams) getDBs() (dbs []DB, err error) {
	return getRegisteredTestDBs(ap.dbNames)
}

func Test_sequence_e2e(t *testing.T) {
	dir := t.TempDir()
	cfg := config.New()

	type testFlags struct {
		in       string
		out      string
		backbone string
		enzymes  []string
		filters  []string
		dbNames  []string
	}

	tests := []testFlags{
		{
			path.Join("..", "..", "test", "input", "backbone.fa"),
			path.Join(dir, "backbone.json"),
			"pSB1A3",
			[]string{"PstI"},
			[]string{},
			[]string{testDB.Name},
		},
		{
			path.Join("..", "..", "test", "input", "BBa_K2224001.fa"),
			path.Join(dir, "BBa_K2224001.json"),
			"pSB1A3",
			[]string{"PstI"},
			[]string{},
			[]string{testDB.Name},
		},
	}

	for _, tt := range tests {
		testInput := createFlagsForTesting(
			tt.in,
			tt.out,
			tt.filters,
			tt.enzymes,
			tt.dbNames)
		testAssemblyParams := &mockAssemblyParams{
			*testInput,
		}

		sols := Sequence(testAssemblyParams, 1, cfg)

		if len(sols) < 1 {
			t.Errorf("no solutions for %s", tt.in)
		}

		for _, s := range sols {
			e := validateJunctions(s, cfg)
			if e != nil {
				t.Logf("failed making %s\n", tt.in)
				t.Error(e)
			}
		}
	}
}

func Test_features(t *testing.T) {
	dir := t.TempDir()

	conf := config.New()
	test1 := createFlagsForTesting(
		"p10 promoter, mEGFP, T7 terminator",
		path.Join(dir, "features.json"),
		[]string{"pSB1A3"},
		[]string{"EcoRI"},
		[]string{testDB.Name},
	)

	type args struct {
		flags *assemblyParamsImpl
		conf  *config.Config
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"test end to end features creation",
			args{
				flags: test1,
				conf:  conf,
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			testAssemblyParams := &mockAssemblyParams{
				assemblyParamsImpl{
					in:           tt.args.flags.in,
					out:          tt.args.flags.out,
					dbNames:      []string{testDB.Name},
					backboneName: tt.args.flags.backboneName,
					enzymeNames:  tt.args.flags.enzymeNames,
					filters:      tt.args.flags.filters,
					identity:     tt.args.flags.identity,
				},
			}

			sols := Features(testAssemblyParams, 1, tt.args.conf)

			if len(sols) < 1 {
				t.Failed()
			}

			for _, s := range sols {
				e := validateJunctions(s, conf)
				if e != nil {
					t.Error(e)
				}
			}
		})
	}
}

func Test_fragments(t *testing.T) {
	c := config.New()
	c.PcrMinFragLength = 10
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
func Test_plasmid_single_plasmid(t *testing.T) {
	dir := t.TempDir()

	c := config.New()

	fs := createFlagsForTesting(
		path.Join("..", "..", "test", "input", "109049.addgene.fa"),
		path.Join(dir, "single-plasmid.json"),
		[]string{},
		[]string{},
		[]string{testDB.Name},
	)

	testAssemblyParams := &mockAssemblyParams{
		*fs,
	}

	assemblies := Sequence(testAssemblyParams, 1, c)

	if !strings.Contains(assemblies[0][0].ID, "109049") {
		t.Fatal("failed to use 109049 to build the plasmid")
	}
}

func createFlagsForTesting(
	in, out string,
	filters, enzymes, dbNames []string,
) *assemblyParamsImpl {

	if strings.Contains(in, ",") {
		inComps := strings.Split(in, ",")
		trimmedComps := []string{}
		for _, entry := range inComps {
			trimmedComps = append(trimmedComps, strings.TrimSpace(entry))
		}
		in = strings.Join(trimmedComps, ",")
	}

	return &assemblyParamsImpl{
		in:          in,
		out:         out,
		dbNames:     dbNames,
		filters:     filters,
		enzymeNames: enzymes,
		identity:    98,
	}
}
