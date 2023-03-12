package repp

import (
	"path"
	"testing"
)

// Test reading of a FASTA file
func Test_read(t *testing.T) {
	type fileRead struct {
		name         string
		file         string
		fragCount    int
		readFeatures bool
	}

	files := []fileRead{
		{
			"113726(circular)",
			path.Join("..", "..", "test", "input", "113726(circular).parent"),
			1,
			false,
		},
		{
			"multi.fasta",
			path.Join("..", "..", "test", "input", "multi.fasta"),
			5,
			false,
		},
		{
			"genbank sequence",
			path.Join("..", "..", "test", "input", "genbank.gb"),
			1,
			false,
		},
		{
			"genbank features",
			path.Join("..", "..", "test", "input", "genbank.gb"),
			66,
			true,
		},
	}

	for _, f := range files {
		fragments, err := read(f.file, f.readFeatures)

		if err != nil {
			t.Error(err)
		}

		if len(fragments) != f.fragCount {
			t.Errorf("failed to load fragments, len=%d, expected=%d", len(fragments), f.fragCount)
		}

		for _, f := range fragments {
			// ensure we got an ID
			if len(f.ID) < 1 {
				t.Error("failed to load an ID for a Frag from FASTA")
			}

			// ensure we got a Seq
			if len(f.Seq) < 1 {
				t.Errorf("failed to parse a sequence for Frag %s", f.ID)
			}
		}
	}
}
