package repp

import (
	"reflect"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

func TestNewFeatureDB(t *testing.T) {
	db := NewFeatureDB()

	if len(db.contents) < 1 {
		t.Fail()
	}
}

func Test_queryFeatures(t *testing.T) {
	tests := []struct {
		name string
		args AssemblyParams
		want [][]string
	}{
		{
			"gather SV40 origin, p10 promoter, mEGFP",
			&assemblyParamsImpl{
				in: "SV40 origin,p10 promoter,mEGFP",
			},
			[][]string{
				{"SV40 origin", "ATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCC"},
				{"p10 promoter", "GACCTTTAATTCAACCCAACACAATATATTATAGTTAAATAAGAATTATTATCAAATCATTTGTATATTAATTAAAATACTATACTGTAAATTACATTTTATTTACAATC"},
				{"mEGFP", "AGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGCGCGGCGAGGGCGAGGGCGATGCCACCAACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTCCTTCAAGGACGACGGCACCTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTTCAACAGCCACAACGTCTATATCACGGCCGACAAGCAGAAGAACGGCATCAAGGCGAACTTCAAGATCCGCCACAACGTCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCAAGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAG"},
			},
		},
		{
			"gather SV40 origin, p10 promoter, mEGFP:rev",
			&assemblyParamsImpl{
				in: "SV40 origin,p10 promoter,mEGFP:rev",
			},
			[][]string{
				{"SV40 origin", "ATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCC"},
				{"p10 promoter", "GACCTTTAATTCAACCCAACACAATATATTATAGTTAAATAAGAATTATTATCAAATCATTTGTATATTAATTAAAATACTATACTGTAAATTACATTTTATTTACAATC"},
				{"mEGFP:REV", "CTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGCTTGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGACGTTGTGGCGGATCTTGAAGTTCGCCTTGATGCCGTTCTTCTGCTTGTCGGCCGTGATATAGACGTTGTGGCTGTTGAAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGGTGCCGTCGTCCTTGAAGGAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTTGGTGGCATCGCCCTCGCCCTCGCCGCGCACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCT"},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			enzymes, err := tt.args.getEnzymes()
			if err != nil {
				t.Fail()
			}
			dbs, err := tt.args.getDBs()
			if err != nil {
				t.Fail()
			}
			backbone, _, err := prepareBackbone(
				tt.args.GetBackboneName(),
				enzymes,
				dbs,
			)
			if err != nil {
				t.Fail()
			}
			if got, _ := queryFeatures(tt.args.GetIn(), backbone, dbs); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("queryFeatures() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_blastFeatures(t *testing.T) {
	type args struct {
		flags          AssemblyParams
		targetFeatures [][]string
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"blast a feature against the part databases",
			args{
				flags: &assemblyParamsImpl{
					identity: 96.0,
				},
				targetFeatures: [][]string{
					{"SV40 origin", "ATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCC"},
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dbs, err := tt.args.flags.getDBs()
			if err != nil {
				t.Fail()
			}
			got := blastFeatures(
				tt.args.flags.GetFilters(),
				tt.args.flags.GetIdentity(),
				dbs,
				tt.args.targetFeatures,
				config.New())

			matches := []match{}
			for _, ms := range got {
				for _, m := range ms {
					matches = append(matches, m.match)
				}
			}

			// confirm that the returned fragments sequences contain at least the full queried sequence
			for _, m := range matches {
				containsTargetSeq := false
				for _, wantedSeq := range tt.args.targetFeatures {
					if ld(m.seq, wantedSeq[1], true) < 5 {
						containsTargetSeq = true
					}
				}

				if !containsTargetSeq {
					t.Fatalf("match with seq %s doesn't contain any of the target features", m.seq)
				}
			}
		})
	}
}
