package repp

import (
	"encoding/json"
	"os"
)

// kv is a simple JSON/serialized key-value store.
// used for enzymes and features right now.
type kv struct {
	contents map[string]string
	path     string
}

func newKV(path string) *kv {
	dat, err := os.ReadFile(path)
	if err != nil {
		rlog.Fatal(err)
	}

	contents := make(map[string]string)
	if err = json.Unmarshal(dat, &contents); err != nil {
		rlog.Fatal(err)
	}

	return &kv{
		contents: contents,
		path:     path,
	}
}

func (k *kv) save() error {
	dat, err := json.MarshalIndent(k.contents, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(k.path, dat, 0644)
}
