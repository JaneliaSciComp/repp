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
