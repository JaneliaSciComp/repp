package repp

import "testing"

func Test_annotate(t *testing.T) {
	type args struct {
		name     string
		seq      string
		output   string
		identity int
		ungapped bool
		dbs      []DB
		filters  []string
		enclosed bool
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"blast a plasmid sequence",
			args{
				"BBa_E0610",
				"TTTACGGCTAGCTCAGTCCTAGGTACAATGCTAGCTACTAGATGAAGTACCTGCTGCCGACCGCGGCGGCGGGTCTGCTGCTGCTGGCGGCGCAGCCGGCGATGGCGGACGATGACGATGACATGAACTTCCCGCGTGCGAGCCGTCTGATGCAGGCGGCGGTGCTGGGTGGCCTGATGGCGGTTAGCGCGGCGGCGACCGCGCAAACCAACCCGTATGCGCGTGGTCCGAACCCGACCGCGGCGAGCCTGGAGGCGAGCGCGGGTCCGTTCACCGTGCGTAGCTTTACCGTTAGCCGTCCGAGCGGTTACGGTGCGGGTACCGTGTACTATCCGACCAACGCGGGTGGCACCGTGGGTGCGATCGCGATTGTTCCGGGTTATACCGCGCGTCAGAGCAGCATCAAATGGTGGGGTCCGCGTCTGGCGAGCCACGGTTTTGTGGTTATCACCATTGATACCAACAGCACCCTGGACCAGCCGAGCAGCCGTAGCAGCCAGCAAATGGCGGCGCTGCGTCAAGTTGCGAGCCTGAACGGTACCAGCAGCAGCCCGATCTACGGCAAGGTGGATACCGCGCGTATGGGCGTTATGGGTTGGAGCATGGGTGGCGGTGGCAGCCTGATTAGCGCGGCGAACAACCCGAGCCTGAAAGCTGCGGCGCCGCAAGCGCCGTGGGACAGCAGCACCAACTTCAGCAGCGTGACCGTTCCGACCCTGATCTTTGCGTGCGAGAACGATAGCATTGCGCCGGTGAACAGCAGCGCGCTGCCGATCTACGACAGCATGAGCCGTAACGCGAAGCAGTTCCTGGAAATTAACGGTGGCAGCCACAGCTGCGCGAACAGCGGTAACAGCAACCAAGCGCTGATTGGCAAGAAAGGTGTGGCGTGGATGAAACGTTTCATGGATAACGACACCCGTTATAGCACCTTTGCGTGCGAAAACCCGAACAGCACCCGTGTTAGCGATTTTCGTACCGCGAATTGCAGCTAATAATACTAGAGAAAGAGGAGAAATACTAGATGAGTGTGATCGCTAAACAAATGACCTACAAGGTTTATATGTCAGGCACGGTCAATGGACACTACTTTGAGGTCGAAGGCGATGGAAAAGGTAAGCCCTACGAGGGGGAGCAGACGGTAAAGCTCACTGTCACCAAGGGCGGACCTCTGCCATTTGCTTGGGATATTTTATCACCACAGTGTCAGTACGGAAGCATACCATTCACCAAGTACCCTGAAGACATCCCTGACTATGTAAAGCAGTCATTCCCGGAGGGCTATACATGGGAGAGGATCATGAACTTTGAAGATGGTGCAGTGTGTACTGTCAGCAATGATTCCAGCATCCAAGGCAACTGTTTCATCTACCATGTCAAGTTCTCTGGTTTGAACTTTCCTCCCAATGGACCTGTCATGCAGAAGAAGACACAGGGCTGGGAACCCAACACTGAGCGTCTCTTTGCACGAGATGGAATGCTGCTAGGAAACAACTTTATGGCTCTGAAGTTAGAAGGAGGCGGTCACTATTTGTGTGAATTTAAAACTACTTACAAGGCAAAGAAGCCTGTGAAGATGCCAGGGTATCACTATGTTGACCGCAAACTGGATGTAACCAATCACAACAAGGATTACACTTCGGTTGAGCAGTGTGAAATTTCCATTGCACGCAAACCTGTGGTCGCCTAATAATACTAGAGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATACGCGGCCGCTTCTAGAGTACTAGTAGCGGCCGCTGCAGTCCGGCAAAAAAGGGCAAGGTGTCACCACCCTGCCCTTTTTCTTTAAAACCGAAAAGATTACTTCGCGTTATGCAGGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCACAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGCTCGAGGCTTGGATTCTCACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCAGATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCCAAGCGAGCTCGATATCAAATTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAGACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGAAATTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCCGTAATATCCAGCTGAACGGTCTGGTTATAGGTACATTGAGCAACTGACTGAAATGCCTCAAAATGTTCTTTACGATGCCATTGGGATATATCAACGGTGGTATATCCAGTGATTTTTTTCTCCATTTTAGCTTCCTTAGCTCCTGAAAATCTCGATAACTCAAAAAATACGCCCGGTAGTGATCTTATTTCATTATGGTGAAAGTTGGAACCTCTTACGTGCCCGATCAACTCGAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCAGAATTTCAGATAAAAAAAATCCTTAGCTTTCGCTAAGGATGATTTCTGG",
				"",
				100,
				false,
				[]DB{testDB},
				[]string{},
				true,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			annotate(tt.args.name, tt.args.seq, tt.args.output, tt.args.identity, tt.args.ungapped, tt.args.dbs, tt.args.filters, tt.args.enclosed, false)
		})
	}
}
