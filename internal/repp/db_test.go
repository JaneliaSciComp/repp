package repp

import (
	"fmt"
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
		Name: "test-db",
		Path: testDbPath,
		Cost: 10,
	}
	testDBs = map[string]DB{
		testDB.Name: testDB,
	}
)

func getRegisteredTestDBs(dbNames []string) (dbs []DB, err error) {
	if len(dbNames) == 0 {
		// if no database was specified - get them all from the manifest
		for _, db := range testDBs {
			dbs = append(dbs, db)
		}
	}

	// filter for matching databases,
	// but only warn the user if a db is not found
	for _, dbName := range dbNames {
		db, ok := testDBs[dbName]
		if ok {
			dbs = append(dbs, db)
		} else {
			rlog.Warnf("DB %s not registered", dbName)
		}
	}

	if len(dbs) == 0 {
		err = fmt.Errorf("None of the requested databases was found.\n The know DBs are: %v", testDBs)
	}

	return
}

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
						Name: "fake_db",
						Path: "/tmp/fake_db.fa",
						Cost: 10,
					},
					{
						Name: "really_fake_db",
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
