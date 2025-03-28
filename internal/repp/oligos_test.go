package repp

import (
	"encoding/csv"
	"reflect"
	"strings"
	"testing"
)

func Test_readOligosFromCSV(t *testing.T) {

	type args struct {
		csvData string
	}
	tests := []struct {
		name      string
		args      args
		want      map[string]oligo
		nextIndex uint
	}{
		{
			name: "parse oligos manifest with header",
			args: args{
				`#comment
				primer_id, sequence
				os1, act
				os2,tgacg`,
			},
			want: map[string]oligo{
				"ACT":   {id: "os1", seq: "act"},
				"TGACG": {id: "os2", seq: "tgacg"},
			},
			nextIndex: 3,
		},
		{
			name: "parse oligos manifest without header",
			args: args{
				`os1,act
				os2,tgacg`,
			},
			want: map[string]oligo{
				"ACT":   {id: "os1", seq: "act"},
				"TGACG": {id: "os2", seq: "tgacg"},
			},
			nextIndex: 3,
		},
		{
			name: "oligo ID has no numeric component",
			args: args{
				`os,act
				os,tgacg
				os,tgacggg`,
			},
			want: map[string]oligo{
				"ACT":     {id: "os", seq: "act"},
				"TGACG":   {id: "os", seq: "tgacg"},
				"TGACGGG": {id: "os", seq: "tgacggg"},
			},
			nextIndex: 4,
		},
		{
			name: "oligo ID has no numeric component and has different prefix",
			args: args{
				`other,act
				other,tgacg
				other,tgacggg`,
			},
			want: map[string]oligo{
				"ACT":     {id: "other", seq: "act"},
				"TGACG":   {id: "other", seq: "tgacg"},
				"TGACGGG": {id: "other", seq: "tgacggg"},
			},
			nextIndex: 1,
		},
		{
			name: "jump in next index",
			args: args{
				`os1,act
				os10,tgacg`,
			},
			want: map[string]oligo{
				"ACT":   {id: "os1", seq: "act"},
				"TGACG": {id: "os10", seq: "tgacg"},
			},
			nextIndex: 11,
		},
		{
			name: "mixed prefixes and unsorted index",
			args: args{
				`os4,act
				os2,tgacg
				other10,actg`,
			},
			want: map[string]oligo{
				"ACT":   {id: "os4", seq: "act"},
				"TGACG": {id: "os2", seq: "tgacg"},
				"ACTG":  {id: "other10", seq: "actg"},
			},
			nextIndex: 5,
		},
		{
			name: "different number of fields per entry but the ones that matter are present",
			args: args{
				`os1,act,,,,,
				os10,tgacg`,
			},
			want: map[string]oligo{
				"ACT":   {id: "os1", seq: "act"},
				"TGACG": {id: "os10", seq: "tgacg"},
			},
			nextIndex: 11,
		},
		{
			name: "row with too few columns",
			args: args{
				`os1,act,,,,,
				os10,
				os11`,
			},
			want: map[string]oligo{
				"ACT": {id: "os1", seq: "act"},
			},
			nextIndex: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := csv.NewReader(strings.NewReader(tt.args.csvData))
			oligos := newOligosDB("oS", false)
			err := readOligosFromCSV(r, oligos)
			if err != nil {
				t.Errorf("%s: Error parsing oligos %v\n", tt.name, err)
			}
			if !reflect.DeepEqual(oligos.indexedOligos, tt.want) {
				t.Errorf("%s: got = %v, want %v\n", tt.name, oligos.indexedOligos, tt.want)
			}
			if oligos.nextOligoID != tt.nextIndex {
				t.Errorf("%s: next index = %d, wanted = %d\n", tt.name, oligos.nextOligoID, tt.nextIndex)
			}
		})
	}
}
