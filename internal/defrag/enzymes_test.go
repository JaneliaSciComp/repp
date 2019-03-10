package defrag

import (
	"reflect"
	"testing"
)

func Test_recogRegex(t *testing.T) {
	type args struct {
		recog string
	}
	tests := []struct {
		name        string
		args        args
		wantDecoded string
	}{
		{
			"decode PpuMI: RGGWCCY",
			args{
				recog: "RGGWCCY",
			},
			"(A|G)GG(A|T)CC(C|T)",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotDecoded := recogRegex(tt.args.recog); gotDecoded != tt.wantDecoded {
				t.Errorf("recogRegex() = %v, want %v", gotDecoded, tt.wantDecoded)
			}
		})
	}

	// should be able to decode every recognition site without failing
	for _, enz := range NewEnzymeDB().enzymes {
		recogRegex(newEnzyme(enz).recog)
	}
}

func Test_digest(t *testing.T) {
	type args struct {
		frag *Frag
		enz  enzyme
	}
	tests := []struct {
		name         string
		args         args
		wantDigested *Frag
		wantBackbone *Backbone
		wantErr      bool
	}{
		{
			"fail with no recognition sequence",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", compCutIndex: 5, seqCutIndex: 1},
			},
			&Frag{},
			&Backbone{},
			true,
		},
		{
			"digest in sequence no overhang",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{name: "TEST", recog: "GAATTC", compCutIndex: 3, seqCutIndex: 3},
			},
			&Frag{
				Seq: "TTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTGAA",
			},
			&Backbone{
				Seq:              "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				Enzyme:           "TEST",
				RecognitionIndex: 24,
				Forward:          true,
			},
			false,
		},
		{
			"digest in reverse complement sequence no overhang",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "CCCAGC", compCutIndex: 3, seqCutIndex: 3}, // rev comp GCTGGG
			},
			&Frag{
				Seq: "GGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGCT",
			},
			&Backbone{
				Seq:              "ATGAGGTTAGCCAAAAAAGCACGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				RecognitionIndex: 22,
				Forward:          false,
			},
			false,
		},
		{
			"digest in sequence positive overhang",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", compCutIndex: 1, seqCutIndex: 5},
			},
			&Frag{
				Seq: "CGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTGAATT",
			},
			&Backbone{
				Seq:              "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				RecognitionIndex: 24,
				Forward:          true,
			},
			false,
		},
		{
			"digest in sequence, negative overhang",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", compCutIndex: 5, seqCutIndex: 1},
			},
			&Frag{
				Seq: "CGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTG",
			},
			&Backbone{
				Seq:              "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				RecognitionIndex: 24,
				Forward:          true,
			},
			false,
		},
		{
			"digest in reverse complement sequence, positive overhang",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "CCCAGC", compCutIndex: 1, seqCutIndex: 5}, // rev comp = GCTGGG
			},
			&Frag{
				Seq: "GGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTGCTGG",
			},
			&Backbone{
				Seq:              "ATGAGGTTAGCCAAAAAAGCACGTGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				RecognitionIndex: 24,
				Forward:          false,
			},
			false,
		},
		{
			"digest in reverse complement sequence, negative overhang",
			args{
				&Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "CCCAGC", compCutIndex: 5, seqCutIndex: 1}, // rev comp = GCTGGG
			},
			&Frag{
				Seq: "GGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTG",
			},
			&Backbone{
				Seq:              "ATGAGGTTAGCCAAAAAAGCACGTGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				RecognitionIndex: 24,
				Forward:          false,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotDigested, gotBackbone, err := digest(tt.args.frag, tt.args.enz)
			if (err != nil) != tt.wantErr {
				t.Errorf("digest() error = %v, wantErr %v", err, tt.wantErr)
				return
			}

			if !reflect.DeepEqual(gotDigested, tt.wantDigested) {
				t.Errorf("digest() = %v, want %v", gotDigested, tt.wantDigested)
			}

			if !reflect.DeepEqual(gotBackbone, tt.wantBackbone) {
				t.Errorf("digest() = %v, want %v", gotBackbone, tt.wantBackbone)
			}
		})
	}
}
