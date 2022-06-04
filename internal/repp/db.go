package repp

import (
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"path"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/spf13/cobra"
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
		names = append(names, db.GetName())
	}
	return names
}

// DB is a single sequence database. It holds the names, sequences, and
// cost of a single sequence source.
type DB struct {
	// Path to the local database in FASTA format.
	Path string `json:"path"`

	// Cost per order from this sequence provider.
	// Eg $65 to order from Addgene.
	Cost float64 `json:"cost"`
}

// GetName returns the name of the database. Just base without an extension right now.
func (d *DB) GetName() string {
	return strings.Replace(path.Base(d.Path), path.Ext(d.Path), "", 1)
}

// AddCmd imports a new sequence database to the REPP directory.
func AddCmd(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("expecting two args: a sequence database and the cost per sequence procurement")
	}

	src := args[0]
	cost := args[1]
	costFloat, err := strconv.ParseFloat(cost, 64)
	if err != nil {
		log.Fatal(err)
	}

	m, err := newManifest()
	if err != nil {
		log.Fatal(err)
	}

	if err = m.add(src, costFloat); err != nil {
		log.Fatal(err)
	}
}

// ListCmd lists the sequence databases and their costs.
func ListCmd(cmd *cobra.Command, args []string) {
	if len(args) > 0 {
		if helperr := cmd.Help(); helperr != nil {
			log.Fatal(helperr)
		}
		log.Fatal("not expecting any arguments")
	}

	m, err := newManifest()
	if err != nil {
		log.Fatal(err)
	}

	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)
	fmt.Fprintf(w, "db\tcost\n")
	for _, db := range m.DBs {
		fmt.Fprintf(w, "%s\t%.2f\n", path.Base(db.Path), db.Cost)
	}
	w.Flush()
}

// DeleteCmd deletes an existing sequence database from the REPP directory.
func DeleteCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		if helperr := cmd.Help(); helperr != nil {
			stderr.Fatal(helperr)
		}
		log.Fatal("expecting one arg: a sequence database name")
	}

	db := args[0]
	m, err := newManifest()
	if err != nil {
		log.Fatal(err)
	}

	if err = m.remove(db); err != nil {
		log.Fatal(err)
	}
}

// newManifest returns a new deserialize Manifest.
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
func (m *manifest) add(from string, cost float64) error {
	src, err := os.Open(from)
	if err != nil {
		return err
	}
	defer src.Close()

	db := DB{
		Path: path.Join(config.SeqDatabaseDir, path.Base(from)),
		Cost: cost,
	}

	dst, err := os.Create(db.Path)
	if err != nil {
		return err
	}
	defer dst.Close()

	_, err = io.Copy(dst, src)
	if err != nil {
		return err
	}

	if err = makeblastdb(db.Path); err != nil {
		return err
	}

	m.DBs[db.GetName()] = db
	return m.save()
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

func dbNames(dbs []DB) (names []string) {
	for _, d := range dbs {
		names = append(names, d.GetName())
	}
	return
}
