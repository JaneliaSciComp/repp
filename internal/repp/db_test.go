package repp

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"
)

var (
	// see test/blast/README.md for a description of where the subfragments
	// in this test fragment's sequence came from (pieces from the 5 fragments)
	// that make up the mock BLAST db
	testDbPath, _ = filepath.Abs(path.Join("..", "..", "test", "db", "db"))
	testDB        = DB{
		Path: testDbPath,
		Cost: 10,
	}
)

func Test_dbNames(t *testing.T) {
	type args struct {
		dbs []DB
	}
	tests := []struct {
		name      string
		args      args
		wantNames []string
	}{
		{
			name: "get db names",
			args: args{
				dbs: []DB{
					{
						Path: "/tmp/fake_db.fa",
						Cost: 10,
					},
					{
						Path: "/tmp/really_fake_db.fa",
						Cost: 15,
					},
				},
			},
			wantNames: []string{"fake_db", "really_fake_db"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotNames := dbNames(tt.args.dbs); !reflect.DeepEqual(gotNames, tt.wantNames) {
				t.Errorf("dbNames() = %v, want %v", gotNames, tt.wantNames)
			}
		})
	}
}
