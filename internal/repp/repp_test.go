package repp

import (
	"os"
	"testing"

	"github.com/Lattice-Automation/repp/internal/config"
)

func TestMain(m *testing.M) {
	config.Setup()
	exitVal := m.Run()
	os.Exit(exitVal)
}
