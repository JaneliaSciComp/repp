// Package seqdb is for interacting with local sequence databases.
//
// DBs are loaded from FASTA files and stored locally along with the cost
// per sequence from the database.
package seqdb

import (
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path"
	"strconv"
	"text/tabwriter"

	"github.com/Lattice-Automation/repp/internal/config"
	"github.com/spf13/cobra"
)

// Manifest is a serializable list of sequence databases.
//
// If sequence databases did not also include meta about cost, this could
// be removed in favor of a simple embedded list of FASTA files/sources.
// This is kept here because of some hypothetical future where the cost
// of sequences depends on their length, shipping method, etc.
type Manifest struct {
	DBs []DB `json:"dbs"`
	path string
}

// DB is a single sequence database. It holds the names, sequences, and
// cost of a single sequence source.
type DB struct {
	// path to the local database in FASTA format.
	Path string `json:"path"`

	// cost per order from this sequence provider.
	// Eg $65 to order from Addgene.
	Cost float64 `json:"cost"`
}

// AddCmd imports a new sequence database to the REPP directory.
func AddCmd(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		cmd.Help()
		log.Fatal("expecting two args: a sequence database and the cost per sequence procurement")
	}

	src := args[0]
	cost := args[1]
	costFloat, err := strconv.ParseFloat(cost, 64)
	if err != nil {
		log.Fatal(err)
	}

	m, err := NewManifest()
	if err != nil {
		log.Fatal(err)
	}

	if err = m.Add(src, costFloat); err != nil {
		log.Fatal(err)
	}
}

// ListCmd lists the sequence databases and their costs.
func ListCmd(cmd *cobra.Command, args []string) {
	if len(args) > 0 {
		cmd.Help()
		log.Fatal("not expecting any arguments")
	}

	m, err := NewManifest()
	if err != nil {
		log.Fatal(err)
	}

	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)
	for _, db := range m.DBs {
		w.Write([]byte(path.Base(db.Path)))
		w.Write([]byte(fmt.Sprintf("%f", db.Cost)))
		w.Flush()
	}
}

// DeleteCmd deletes an existing sequence database from the REPP directory.
func DeleteCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		cmd.Help()
		log.Fatal("expecting one arg: a sequence database name")
	}

	db := args[1]
	m, err := NewManifest()
	if err != nil {
		log.Fatal(err)
	}

	if err = m.Remove(db); err != nil {
		log.Fatal(err)
	}
}

// NewManifest returns a new deserialize Manifest.
func NewManifest() (*Manifest, error) {
	contents, err := ioutil.ReadFile(config.SeqDatabaseManifest)
	if err != nil {
		if os.IsNotExist(err) {
			return &Manifest{
				path: config.SeqDatabaseManifest,
				DBs: []DB{},
			}, nil
		}
		return nil, err
	}

	var manifest *Manifest
	if err = json.Unmarshal(contents, manifest); err != nil {
		return nil, err
	}
	return manifest, nil
}

// Add imports a FASTA sequence database into REPP, storing it in the manifest.
func (m *Manifest) Add(from string, cost float64) error {
	src, err := os.Open(from)
	if err != nil {
		return err
	}

	to := path.Join(m.path, path.Base(from))
	dst, err := os.Open(to)
	if err != nil {
		return err
	}
	_, err = io.Copy(src, dst)
	if err != nil {
		return err
	}

	m.DBs = append(m.DBs, DB{
		Path: to,
		Cost: cost,
	})
	return m.save()
}

// Remove deletes a local, repp-managed FASTA file and removes it from the manifest
func (m *Manifest) Remove(name string) error {
	dbs := []DB{}
	for _, db := range m.DBs {
		if path.Base(db.Path) == name {
			if err := os.Remove(db.Path); err != nil {
				return err
			}
		} else {
			dbs = append(dbs, db)
		}
	}
	if len(dbs) == len(m.DBs) {
		return fmt.Errorf("no DB found with name: %s", name)
	}

	m.DBs = dbs
	return m.save()
}

func (m *Manifest) save() error {
	contents, err := json.MarshalIndent(m, "", "  ")
	if err != nil {
		return err
	}
	return ioutil.WriteFile(m.path, contents, 0644)
}
