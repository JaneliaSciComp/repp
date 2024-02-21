package repp

import (
	"fmt"
	"math"
	"reflect"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

func Test_assembly_add(t *testing.T) {
	c := config.New()

	c.FragmentsMaxCount = 5
	c.PcrPrimerMaxEmbedLength = 0
	c.PcrMinFragLength = 0
	c.SyntheticMaxLength = 100
	c.SyntheticFragmentCost = map[int]config.SynthCost{
		100: {
			Fixed: true,
			Cost:  1.0,
		},
	}
	c.SyntheticFragmentFactor = 2

	sl := 100

	n1 := &Frag{
		ID:       "1",
		uniqueID: "1",
		fragType: pcr,
		start:    0,
		end:      50,
		conf:     c,
	}
	n2 := &Frag{
		ID:       "2",
		uniqueID: "2",
		fragType: pcr,
		start:    20,
		end:      80,
		conf:     c,
	}
	n3 := &Frag{
		ID:       "3",
		uniqueID: "3",
		fragType: pcr,
		start:    60,
		end:      100,
		conf:     c,
	}
	altn3 := &Frag{
		ID:       n3.ID,
		uniqueID: n3.uniqueID,
		fragType: pcr,
		start:    n3.start + c.SyntheticMaxLength,
		end:      n3.end + c.SyntheticMaxLength,
		conf:     c,
	}
	n4 := &Frag{
		ID:       "1",
		uniqueID: "1",
		fragType: pcr,
		start:    100,
		end:      150,
		conf:     c,
	}

	// create the frags for testing
	type fields struct {
		frags        []*Frag
		cost         float64
		adjustedCost float64
		synths       int
	}
	type args struct {
		n *Frag
	}

	createAssemblyFrom := func(fs []*Frag,
		selfAnnealing bool,
		costFunc func() (float64, float64),
		synths int) assembly {
		cost, adjustedCost := costFunc()

		return assembly{
			frags:         fs,
			selfAnnealing: selfAnnealing,
			cost:          cost,
			adjustedCost:  adjustedCost,
			synths:        synths,
		}
	}
	tests := []struct {
		name            string
		fields          fields
		args            args
		wantNewAssembly assembly
		wantComplete    bool
		wantErr         error
	}{
		{
			"add with overlap",
			fields{
				frags:  []*Frag{n1},
				cost:   0,
				synths: 0,
			},
			args{
				n: n2,
			},
			createAssemblyFrom([]*Frag{n1, n2},
				false,
				func() (float64, float64) {
					n1c, n1ac := n1.cost(true)
					n1Ton2c, n1Ton2ac := n1.costTo(n2)
					return n1c + n1Ton2c, n1ac + n1Ton2ac
				},
				0),
			false,
			nil,
		},
		{
			"add with synthesis",
			fields{
				frags:        []*Frag{n1},
				cost:         10.0,
				adjustedCost: 10.0,
				synths:       0,
			},
			args{
				n: n3,
			},
			createAssemblyFrom([]*Frag{n1, n3},
				false,
				func() (float64, float64) {
					n1c, n1ac := n1.cost(true)
					c, ac := n1.costTo(n3)
					return 10 + n1c + c, 10 + n1ac + ac
				},
				1),
			true,
			nil,
		},
		{
			"add with completion/circularization",
			fields{
				frags:  []*Frag{n1, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n4,
			},
			createAssemblyFrom([]*Frag{n1, n2, n3},
				true,
				func() (float64, float64) {
					n3c, n3ac := n3.cost(true)
					c, ac := n1.costTo(n3)
					return 10. + n3c + c, n3ac + ac
				},
				0),
			true,
			nil,
		},
		{
			"add with completion requiring synthesis",
			fields{
				frags:        []*Frag{n1, n2},
				cost:         16.4,
				adjustedCost: 18.4,
				synths:       0,
			},
			args{
				// a Frag that's too far away for straightforward annealing
				n: altn3,
			},
			createAssemblyFrom([]*Frag{n1, n2, altn3},
				false,
				func() (float64, float64) {
					n3c, n3ac := altn3.cost(false)
					return n3c + 48.4, n3ac + 52.4
				},
				1),
			true,
			nil,
		},
		{
			"don't exceed fragment limit",
			fields{
				frags:  []*Frag{n1, n2, n3, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n2,
			},
			assembly{},
			false,
			fmt.Errorf("the resulted assembly has  more fragments than allowed (6 > 5)"),
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := assembly{
				frags:        tt.fields.frags,
				cost:         tt.fields.cost,
				adjustedCost: tt.fields.adjustedCost,
				synths:       tt.fields.synths,
			}
			gotNewAssembly, gotComplete, gotErr := extendAssembly(a, tt.args.n, 5, sl, false)
			if !reflect.DeepEqual(gotNewAssembly.frags, tt.wantNewAssembly.frags) ||
				math.Abs(gotNewAssembly.cost-tt.wantNewAssembly.cost) > 0.1 ||
				math.Abs(gotNewAssembly.adjustedCost-tt.wantNewAssembly.adjustedCost) > 0.1 {
				t.Errorf("assembly.add() gotNewAssembly = %v, want %v", gotNewAssembly, tt.wantNewAssembly)
			}
			if gotComplete != tt.wantComplete {
				t.Errorf("assembly.add() gotComplete = %v, want %v", gotComplete, tt.wantComplete)
			}
			if gotErr == nil && tt.wantErr != nil ||
				gotErr != nil && tt.wantErr == nil ||
				gotErr != nil && tt.wantErr != nil && gotErr.Error() != tt.wantErr.Error() {
				t.Errorf("assembly.add() gotComplete = %v, want %v", gotErr, tt.wantErr)
			}
		})
	}
}
func Test_assembly_len(t *testing.T) {
	c := config.New()

	n1 := &Frag{
		uniqueID: "1",
		start:    0,
		end:      50,
		conf:     c,
	}
	n2 := &Frag{
		uniqueID: "2",
		start:    20,
		end:      80,
		conf:     c,
	}

	type fields struct {
		frags  []*Frag
		cost   float64
		synths int
	}
	tests := []struct {
		name   string
		fields fields
		want   int
	}{
		{
			"length without synths",
			fields{
				frags:  []*Frag{n1, n2},
				synths: 0,
			},
			2,
		},
		{
			"length with synths",
			fields{
				frags:  []*Frag{n1, n2},
				synths: 2,
			},
			4,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				frags:  tt.fields.frags,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			if got := a.len(); got != tt.want {
				t.Errorf("assembly.len() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_assembly_duplicates(t *testing.T) {
	type args struct {
		frags       []*Frag
		minHomology int
		maxHomology int
	}
	tests := []struct {
		name         string
		args         args
		want         bool
		wantJunction string
	}{
		{
			"no false positive",
			args{
				frags: []*Frag{
					{
						Seq: "ATACCTACTATGGATGACGTAGCAAC",
					},
					{
						Seq: "AGCAACTCGTTGATATCCACGTA",
					},
					{
						Seq: "CCACGTAGGTGCATGATGAGATGA",
					},
					{
						Seq: "TGAGATGATCTACTGTATACCTACT",
					},
				},
				minHomology: 5,
				maxHomology: 10,
			},
			false,
			"",
		},
		{
			"assembly with a self-annealing Frag",
			args{
				frags: []*Frag{
					{
						Seq: "CAGATGACGATGGCAACTGAGATGAGACCAGATGACGATG", // <- Frag (if much larger) has the chance to circularize
					},
					{
						Seq: "CAGATGACGATGTCGTTGATATACCTACTGGAGAGCACAG",
					},
					{
						Seq: "TGGAGAGCACAGATGGATGACGTAATGATGATGACCGCAAC",
					},
					{
						Seq: "ACCGCAACTCGTTGATATACCTACTCAGATGACGAT",
					},
				},
				minHomology: 5,
				maxHomology: 20,
			},
			true,
			"CAGATGACGATG",
		},
		{
			"assembly with a duplicate junction",
			args{
				frags: []*Frag{
					{
						Seq:   "ATGATGCCACGTGCAACTGAGATGAGACCAGATGACGATG", // <- same junction
						start: 0,
					},
					{
						Seq:   "CAGATGACGATGTCGTTGATATACCTACTGGAGAGCACAG",
						start: 0,
					},
					{
						Seq:   "TGGAGAGCACAGATGGATGACGTAATGACAGATGACGATG", // <- same junction
						start: 0,
					},
					{
						Seq:   "CAGATGACGATGACCGCAACTCGTTGATGATGCCAC",
						start: 0,
					},
				},
				minHomology: 5,
				maxHomology: 20,
			},
			true,
			"CAGATGACGATG",
		},
		{
			"another false positive to avoid",
			args{
				frags: []*Frag{
					{
						Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG", // this shouldn't be flagged as anything
					},
					{
						Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					},
					{
						Seq: "TCGACTAGCTAGAACTGATGCTAGACGTGCTAGCTACA",
					},
				},
				minHomology: 8,
				maxHomology: 20,
			},
			false,
			"",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			isDuplicate, _, _, duplicateSeq := duplicates(tt.args.frags, tt.args.minHomology, tt.args.maxHomology)

			if isDuplicate != tt.want {
				t.Errorf("assembly.duplicates() = %v, want %v", isDuplicate, tt.want)
			}

			if duplicateSeq != tt.wantJunction {
				t.Errorf("assembly.duplicates() = %v, want %v", duplicateSeq, tt.wantJunction)
			}
		})
	}
}
