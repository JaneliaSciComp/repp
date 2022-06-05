package repp

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"path"
	"strconv"
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

// AddCmd imports a new sequence database to the REPP directory.
func AddCmd(cmd *cobra.Command, args []string) {
	_, err := os.Stdin.Stat()
	if err != nil {
		if helperr := cmd.Help(); helperr != nil {
			rlog.Fatal(helperr)
		}
		rlog.Fatal("no stdin passed to 'repp add database'. See example.")
	}

	name := cmd.Flag("name").Value.String()
	cost := cmd.Flag("cost").Value.String()
	costFloat, err := strconv.ParseFloat(cost, 64)
	if err != nil {
		rlog.Fatal(err)
	}

	m, err := newManifest()
	if err != nil {
		rlog.Fatal(err)
	}

	if err = m.add(name, costFloat); err != nil {
		rlog.Fatal(err)
	}
}

// ListCmd lists the sequence databases and their costs.
func ListCmd(cmd *cobra.Command, args []string) {
	m, err := newManifest()
	if err != nil {
		rlog.Fatal(err)
	}

	if m.empty() {
		rlog.Fatal("No databases loaded. See 'repp add database'")
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
	db := args[0]
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
func (m *manifest) add(name string, cost float64) error {
	db := DB{
		Name: name,
		Path: path.Join(config.SeqDatabaseDir, name),
		Cost: cost,
	}

	l := rlog.With("path", db.Path, "name", name, "cost", cost)
	dst, err := os.Create(db.Path)
	if err != nil {
		l.Error("failed to create db")
		return err
	}
	l.Debug("created db")
	defer dst.Close()

	reader := bufio.NewReader(os.Stdin)
	n, err := reader.WriteTo(dst)
	if err != nil {
		l.Errorw("failed copying stdin", "err", err)
		return err
	}
	l.Debugw("copied stdin", "bytes", n)
	if n == 0 {
		l.Fatal("no stdin to create db")
	}

	if err = makeblastdb(db.Path); err != nil {
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

func dbNames(dbs []DB) (names []string) {
	for _, d := range dbs {
		names = append(names, d.Name)
	}
	return
}
