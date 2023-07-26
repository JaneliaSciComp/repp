package repp

import (
	"reflect"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

func Test_assembly_add(t *testing.T) {
	c := config.New()

	c.FragmentsMaxCount = 5
	c.PcrPrimerMaxEmbedLength = 0
	c.PcrMinLength = 0
	c.SyntheticMaxLength = 100
	c.SyntheticFragmentCost = map[int]config.SynthCost{
		100: {
			Fixed: true,
			Cost:  1.0,
		},
	}
	c.SyntheticFragmentPenalty = 2

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
		costFunc func() (float64, float64),
		synths int) assembly {
		cost, adjustedCost := costFunc()

		return assembly{
			frags:        fs,
			cost:         cost,
			adjustedCost: adjustedCost,
			synths:       synths,
		}
	}
	tests := []struct {
		name            string
		fields          fields
		args            args
		wantNewAssembly assembly
		wantCreated     bool
		wantComplete    bool
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
				func() (float64, float64) {
					return n1.costTo(n2)
				},
				0),
			true,
			false,
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
				func() (float64, float64) {
					c, ac := n1.costTo(n3)
					return 10 + c, 10 + ac
				},
				1),
			true,
			true,
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
				func() (float64, float64) {
					return 10, 0
				},
				0),
			true,
			true,
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
				n: &Frag{
					ID:       n3.ID,
					uniqueID: n3.uniqueID,
					fragType: pcr,
					start:    n3.start + c.SyntheticMaxLength,
					end:      n3.end + c.SyntheticMaxLength,
					conf:     c,
				},
			},
			createAssemblyFrom([]*Frag{n1, n2, {
				ID:       n3.ID,
				uniqueID: n3.uniqueID,
				fragType: pcr,
				start:    n3.start + c.SyntheticMaxLength,
				end:      n3.end + c.SyntheticMaxLength,
				conf:     c,
			}},
				func() (float64, float64) {
					return 48.4, 52.4
				},
				1),
			true,
			true,
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
			false,
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
			gotNewAssembly, gotCreated, gotComplete := createNewAssembly(a, tt.args.n, 5, sl, false)
			if !reflect.DeepEqual(gotNewAssembly, tt.wantNewAssembly) {
				t.Errorf("assembly.add() gotNewAssembly = %v, want %v", gotNewAssembly, tt.wantNewAssembly)
			}
			if gotCreated != tt.wantCreated {
				t.Errorf("assembly.add() gotCreated = %v, want %v", gotCreated, tt.wantCreated)
			}
			if gotComplete != tt.wantComplete {
				t.Errorf("assembly.add() gotComplete = %v, want %v", gotComplete, tt.wantComplete)
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

func Test_countMaps(t *testing.T) {
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
	n3 := &Frag{
		uniqueID: "3",
		start:    60,
		end:      100,
		conf:     c,
	}

	a1 := assembly{
		frags: []*Frag{
			n1, n1,
		},
		cost: 11.0,
	}
	a2 := assembly{
		frags: []*Frag{
			n1, n2, n1,
		},
		cost: 12.5,
	}
	a3 := assembly{
		frags: []*Frag{
			n2, n3, n2,
		},
		cost: 12.0,
	}
	a4 := assembly{
		frags: []*Frag{
			n1, n2, n3, n1,
		},
		cost: 10.0,
	}
	a5 := assembly{
		frags: []*Frag{
			n2, n3, n1, n2,
		},
		cost: 10.5,
	}

	type args struct {
		assemblies []assembly
	}
	tests := []struct {
		name          string
		args          args
		wantParetoSet map[int][]assembly
	}{
		{
			"gen pSet up to 3",
			args{
				assemblies: []assembly{a1, a2, a3, a4, a5},
			},
			map[int][]assembly{
				2: {a1},
				3: {a3, a2},
				4: {a4, a5},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if _, gotParetoSet := groupAssembliesByCount(tt.args.assemblies); !reflect.DeepEqual(gotParetoSet, tt.wantParetoSet) {
				t.Errorf("pareto() = %v, want %v", gotParetoSet, tt.wantParetoSet)
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
