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
// be removed in favor of a simple directory of FASTA files/sources.
// This is here because of some hypothetical future where the cost
// of sequences depends on their length, shipping method, etc.
type manifest struct {
	DBs []DB `json:"dbs"`
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
		cmd.Help()
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
		cmd.Help()
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
				DBs: []DB{},
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
	// TODO: quiteee the hack to avoid duplicates
	m.remove(path.Base(from))

	src, err := os.Open(from)
	if err != nil {
		return err
	}
	defer src.Close()

	toBase := strings.Replace(path.Base(from), path.Ext(from), "", 1)
	to := path.Join(config.SeqDatabaseDir, toBase)
	dst, err := os.Create(to)
	if err != nil {
		return err
	}
	defer dst.Close()

	_, err = io.Copy(dst, src)
	if err != nil {
		return err
	}

	if err = makeblastdb(to); err != nil {
		return err
	}

	m.DBs = append(m.DBs, DB{
		Path: to,
		Cost: cost,
	})
	return m.save()
}

// remove deletes a local, repp-managed FASTA file and removes it from the manifest
func (m *manifest) remove(name string) error {
	dbs := []DB{}
	for _, db := range m.DBs {
		if path.Base(db.Path) == name {
			// TODO: also remove the other BLAST files, rather than just deleting the FASTA
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

func (m *manifest) save() error {
	contents, err := json.MarshalIndent(m, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(config.SeqDatabaseManifest, contents, 0644)
}
