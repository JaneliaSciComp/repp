package repp

import (
	"os"
	"testing"
)

func Test_writeGenbank(t *testing.T) {
	dir := t.TempDir()
	output, err := os.CreateTemp(dir, "*.gb")
	if err != nil {
		t.Error(err)
	}

	type args struct {
		filename string
		name     string
		seq      string
		frags    []*Frag
		feats    []match
	}
	tests := []struct {
		name string
		args args
	}{
		{
			name: "include forward and reverse fields",
			args: args{
				output.Name(),
				"mock part",
				"aattgtgagcggataacaattgacattgtgagcggataacaagatactgagcacatactagagaaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagaggcctgctgcaaacgacgaaaactacgctttagtagcttaataatactagagtcacactggctcaccttcgggtgggcctttctgcgtttatatactagagagagaatataaaaagccagattattaatccggcttttttattattt",
				[]*Frag{},
				[]match{
					{
						entry:               "feature 1",
						queryStart:          0,
						queryEnd:            10,
						queryRevCompMatch:   false,
						subjectRevCompMatch: false,
					},
					{
						entry:               "feature 2",
						queryStart:          15,
						queryEnd:            20,
						queryRevCompMatch:   false,
						subjectRevCompMatch: true,
					},
				},
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			writeGenbank(tt.args.filename, tt.args.name, tt.args.seq, tt.args.frags, tt.args.feats)
		})
	}
}
