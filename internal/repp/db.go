package repp

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"path"
	"text/tabwriter"

	"github.com/Lattice-Automation/repp/internal/config"
)

// manifest is a serializable list of sequence databases.
//
// If sequence databases did not also include meta about cost, this could
// be removed in favor of a simple directory of FASTA files (1 per database).
//
// This is here because of some hypothetical future where the cost
// of sequences depends on their length, shipping method, etc.
type manifest struct {
	// DBs is a map from DB name (base of originally added DB file) to DB
	DBs map[string]DB `json:"dbs"`
}

// GetNames returns the list of known DB names.
func (m *manifest) GetNames() (names []string) {
	for _, db := range m.DBs {
		names = append(names, db.Name)
	}
	return names
}

// DB is a single sequence database. It holds the names, sequences, and
// cost of a single sequence source.
type DB struct {
	// Name of the db
	Name string `json:"name"`

	// Path to the local database in FASTA format.
	Path string `json:"path"`

	// Cost per order from this sequence provider.
	// Eg $65 to order from Addgene.
	Cost float64 `json:"cost"`
}

// AddDatabase imports one or more sequence files into a BLAST database to the REPP directory.
func AddDatabase(dbName string, seqFiles []string, cost float64, appendMode bool) (err error) {
	dbSequenceFilepath := path.Join(config.SeqDatabaseDir, dbName)

	var openFlags int
	if appendMode {
		openFlags = os.O_WRONLY | os.O_CREATE | os.O_APPEND
	} else {
		openFlags = os.O_WRONLY | os.O_CREATE
	}
	dbSeqFile, err := os.OpenFile(dbSequenceFilepath, openFlags, 0644)
	if err != nil {
		rlog.Errorf("Error creating database sequence file %f\n", dbSequenceFilepath)
		return err
	}
	defer dbSeqFile.Close()

	if len(seqFiles) == 0 {
		// try to read from stdin
		_, err := os.Stdin.Stat()
		if err != nil {
			rlog.Warnf("Error reading sequence from the standard input")
			return err
		}
		dbSeqInput := os.Stdin
		dbSeqReader := bufio.NewReader(dbSeqInput)

		if _, err = dbSeqReader.WriteTo(dbSeqFile); err != nil {
			rlog.Errorf("Error writing database sequence to %f\n", dbSequenceFilepath)
			return err
		}
	} else {
		dbSeqs, err := multiFileRead(seqFiles)
		if err != nil {
			rlog.Errorf("Error reading one or more sequence files into the database")
			return err
		}
		err = writeFragsToFastaFile(dbSeqs, dbSeqFile)
		if err != nil {
			rlog.Errorf("Error writing database sequence to %f\n", dbSequenceFilepath)
			return err
		}
	}

	m, err := newManifest()
	if err != nil {
		rlog.Fatal(err)
	}

	if err = m.add(dbName, dbSequenceFilepath, cost); err != nil {
		rlog.Fatal(err)
	}

	return
}

// ListCmd lists the sequence databases and their costs.
func ListDatabases() {
	m, err := newManifest()
	if err != nil {
		rlog.Fatal(err)
	}

	if m.empty() {
		rlog.Fatal("No databases loaded. See 'repp add database'")
	}

	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)
	fmt.Fprintf(w, "name\tcost\n")
	for _, db := range m.DBs {
		fmt.Fprintf(w, "%s\t%.2f\n", path.Base(db.Path), db.Cost)
	}
	w.Flush()
}

// DeleteCmd deletes an existing sequence database from the REPP directory.
func DeleteDatabase(db string) {
	m, err := newManifest()
	if err != nil {
		rlog.Fatal(err)
	}

	if err = m.remove(db); err != nil {
		rlog.Fatal(err)
	}
}

// newManifest returns a new deserialized Manifest.
func newManifest() (*manifest, error) {
	contents, err := os.ReadFile(config.SeqDatabaseManifest)
	if err != nil {
		if os.IsNotExist(err) {
			return &manifest{
				DBs: map[string]DB{},
			}, nil
		}
		return nil, err
	}

	m := &manifest{}
	if err = json.Unmarshal(contents, m); err != nil {
		return nil, err
	}
	return m, nil
}

// add imports a FASTA sequence database into REPP, storing it in the manifest.
func (m *manifest) add(dbName string, seqFilepath string, cost float64) error {
	db := DB{
		Name: dbName,
		Path: seqFilepath,
		Cost: cost,
	}
	l := rlog.With("path", db.Path, "name", dbName, "cost", cost)
	if err := makeblastdb(db.Path); err != nil {
		l.Error("failed to makeblastdb")
		return err
	}
	l.Debug("ran makeblastdb")

	m.DBs[db.Name] = db

	return m.save()
}

// empty returns whether the manifest lacks any database
func (m *manifest) empty() bool {
	return len(m.DBs) == 0
}

// remove deletes a local, repp-managed FASTA file and removes it from the manifest
func (m *manifest) remove(name string) error {
	if _, ok := m.DBs[name]; !ok {
		return fmt.Errorf("no DB found with name %s", name)
	}
	delete(m.DBs, name)
	return m.save()
}

func (m *manifest) save() error {
	contents, err := json.MarshalIndent(m, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(config.SeqDatabaseManifest, contents, 0644)
}

func getRegisteredDBs(dbNames []string) (dbs []DB, err error) {
	m, err := newManifest()
	if err != nil {
		rlog.Fatalf("failed to get DB manifest: %v", err)
	}

	if len(dbNames) == 0 {
		// if no database was specified - get them all from the manifest
		for _, db := range m.DBs {
			dbs = append(dbs, db)
		}
		return
	}

	// filter for matching databases,
	// but only warn the user if a db is not found
	for _, dbName := range dbNames {
		db, ok := m.DBs[dbName]
		if ok {
			dbs = append(dbs, db)
		} else {
			rlog.Warnf("DB %s not registered", dbName)
		}
	}

	if len(dbs) == 0 {
		err = fmt.Errorf("None of the requested databases was found.\n The know DBs are: %v", m.GetNames())
	}

	return
}

func dbNames(dbs []DB) (names []string) {
	for _, d := range dbs {
		names = append(names, d.Name)
	}
	return
}
