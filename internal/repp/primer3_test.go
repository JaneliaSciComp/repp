package repp

import (
	"math"
	"reflect"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

func Test_primer3_shrink(t *testing.T) {
	type args struct {
		seq    string
		config *config.Config
		last   *Frag
		n      *Frag
		next   *Frag
	}
	tests := []struct {
		name string
		args args
		want *Frag
	}{
		{
			"shrink Frag with an excessive amount of homology",
			args{
				seq: "", // not important here
				config: func() *config.Config {
					c := config.New()
					c.FragmentsMaxHomology = 10 // much less than normal
					c.PcrMinFragLength = 20
					c.PcrPrimerUseStrictConstraints = false
					return c
				}(),
				last: &Frag{
					start: 0,
					end:   100,
				},
				n: &Frag{
					Seq:   "GGGGGAACGCTGAAGATCTCTTCTTCTCATGACTGAACTCGCGAGGGTCGTGATGTCGGTTCCTTCAAAGGTTAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTCTTGAATCCTCGGTCCCCCTTGTCTTTCCAGATTAATCCATTTCCCTCATTCACGAGCTTACCAAGTCAACATTGGTATATGAATGCGACCTTGAAGAGGCCGCTTAAAAATGGCAGTGGTTGAT",
					start: 90,
					end:   300,
				},
				next: &Frag{
					start: 250,
					end:   500,
				},
			},
			&Frag{
				Seq:   "GGGGGAACGCTGAAGATCTCTTCTTCTCATGACTGAACTCGCGAGGGTCGTGATGTCGGTTCCTTCAAAGGTTAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTCTTGAATCCTCGGTCCCCCTTGTCTTTCCAGATTAATCCATTTCCCTCATTCACGAGCTTACCAAGTCAACATTGGTATATGAAT",
				start: 90,
				end:   260,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := newPrimer3(tt.args.seq, tt.args.config)
			if got := p.shrink(tt.args.last, tt.args.n, tt.args.next); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("primer3.shrink() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_bpToAdd(t *testing.T) {
	c := config.New()
	c.PcrPrimerMaxEmbedLength = 20
	c.FragmentsMinHomology = 10

	p := primer3{
		config: c,
	}

	type args struct {
		left  *Frag
		right *Frag
	}
	tests := []struct {
		name        string
		args        args
		wantBpToAdd int
	}{
		{
			"no added homology is needed",
			args{
				left: &Frag{
					start: 0,
					end:   20,
					conf:  c,
				},
				right: &Frag{
					start: 10,
					end:   30,
					conf:  c,
				},
			},
			0,
		},
		{
			"no added homology is needed - lots of overlap",
			args{
				left: &Frag{
					start: 0,
					end:   50,
					conf:  c,
				},
				right: &Frag{
					start: 10,
					end:   30,
					conf:  c,
				},
			},
			0,
		},
		{
			"add homology to each Frag",
			args{
				left: &Frag{
					start: 0,
					end:   10,
					conf:  c,
				},
				right: &Frag{
					start: 16,
					end:   30,
					conf:  c,
				},
			},
			12,
		},
		{
			"correct bp to share when negative",
			args{
				left: &Frag{
					start: -50,
					end:   20,
					conf:  c,
				},
				right: &Frag{
					start: 0,
					end:   1050,
					conf:  c,
				},
			},
			0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotBpToAdd := p.bpToAdd(tt.args.left, tt.args.right); gotBpToAdd != tt.wantBpToAdd {
				t.Errorf("bpToAdd() = %v, want %v", gotBpToAdd, tt.wantBpToAdd)
			}
		})
	}
}

func Test_mutatePrimers(t *testing.T) {
	type args struct {
		n        *Frag
		Seq      string
		addLeft  int
		addRight int
	}
	tests := []struct {
		name        string
		args        args
		wantMutated *Frag
	}{
		{
			"add homology to both sides of a Frag",
			args{
				n: &Frag{
					Seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
					start: 10,
					end:   39,
					Primers: []Primer{
						{
							Seq: "TGACCTCGGC",
							Range: ranged{
								start: 10,
								end:   20,
							},
							Strand: true,
						},
						{
							Seq: "CGCCGTAGTA", // rev comp TACTACGGCG
							Range: ranged{
								start: 30,
								end:   39,
							},
							Strand: false,
						},
					},
				},
				Seq:      "GATCACTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTTGGAG",
				addLeft:  5,
				addRight: 6,
			},
			&Frag{
				PCRSeq: "CTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTT",
				Seq:    "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				start:  10,
				end:    39,
				Primers: []Primer{
					{
						Seq: "CTCGATGACCTCGGC",
						Range: ranged{
							start: 5,
							end:   20,
						},
						Strand: true,
					},
					{
						Seq: "AAGAATCGCCGTAGTA", //  rev comp TACTACGGCGATTCTT
						Range: ranged{
							start: 30,
							end:   45,
						},
						Strand: false,
					},
				},
			},
		},
		{
			"add nothing",
			args{
				n: &Frag{
					Seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
					start: 10,
					end:   39,
					Primers: []Primer{
						{
							Seq: "TGACCTCGGC",
							Range: ranged{
								start: 10,
								end:   20,
							},
							Strand: true,
						},
						{
							Seq: "CGCCGTAGTA", // rev comp TACTACGGCG
							Range: ranged{
								start: 30,
								end:   39,
							},
							Strand: false,
						},
					},
				},
				Seq:      "GATCACTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTTGGAG",
				addLeft:  0,
				addRight: 0,
			},
			&Frag{
				Seq:    "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				PCRSeq: "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				start:  10,
				end:    39,
				Primers: []Primer{
					{
						Seq: "TGACCTCGGC",
						Range: ranged{
							start: 10,
							end:   20,
						},
						Strand: true,
					},
					{
						Seq: "CGCCGTAGTA", // rev comp TACTACGGCG
						Range: ranged{
							start: 30,
							end:   39,
						},
						Strand: false,
					},
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotMutated := mutatePrimers(tt.args.n, tt.args.Seq, tt.args.addLeft, tt.args.addRight); !reflect.DeepEqual(gotMutated, tt.wantMutated) {
				t.Errorf("mutateNodePrimers() = %v, want %v", gotMutated, tt.wantMutated)
			}
		})
	}
}

// these estimated hairpin tms jump around when the primer3 version changes
func Test_hairpin(t *testing.T) {
	c := config.New()

	type args struct {
		seq  string
		conf *config.Config
	}
	tests := []struct {
		name     string
		args     args
		wantMelt float64
	}{
		{
			name: "find hairpin of ~85 degrees",
			args: args{
				"TGTGCACTCATCATCATCATCGGGGGGGGGGGGTGAACACTATCCCCCCCCCCCCCCA",
				c,
			},
			wantMelt: 85.0,
		},
		{
			name: "return 0 when no hairpin found",
			args: args{
				"TGTGcactcatcatcCCCA",
				c,
			},
			wantMelt: 0.0,
		},
		{
			name: "return the right-most hairpin when >60bp",
			args: args{
				"TGTGcactcatcatcaacacaactacgtcgatcagctacgatcgatcgatgctgatcgatatttatatcgagctagctacggatcatcGGGGGGGGGGGGTGAACACTATCCCCCCCCCCCCCCA",
				c,
			},
			wantMelt: 85.0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotMelt := hairpin(tt.args.seq, tt.args.conf); math.Abs(gotMelt-tt.wantMelt) > 10 {
				t.Errorf("hairpin() = %v, want %v", gotMelt, tt.wantMelt)
			}
		})
	}
}
