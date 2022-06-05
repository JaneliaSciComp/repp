// Package config is for app wide settings
package config

import (
	"embed"
	"log"
	"math"
	"os"
	"path"
	"path/filepath"
	"sort"

	"github.com/mitchellh/go-homedir"
	"github.com/mitchellh/mapstructure"
	"github.com/spf13/viper"
	"gopkg.in/yaml.v2"
)

var (
	home, _ = homedir.Dir()

	// reppDir is the root directory where repp settings and database files live
	reppDir = filepath.Join(home, ".repp")

	// configPath is the path to a local/default config file
	configPath = filepath.Join(reppDir, "config.yaml")

	// Primer3Config is the path to a primer3 config directory
	// primer3 is (overly) particular about the trailing slash
	Primer3Config = filepath.Join(reppDir, "primer3_config") + string(os.PathSeparator)

	// FeatureDB is the path to the features file
	FeatureDB = filepath.Join(reppDir, "features.json")

	// EnzymeDB is the path to the enzymes file
	EnzymeDB = filepath.Join(reppDir, "enzymes.json")

	// SeqDatabaseDir is the path to a directory of sequence databases.
	SeqDatabaseDir = filepath.Join(reppDir, "dbs")

	// SeqDatabaseManifest is the path to the manifest file for the sequence databases.
	SeqDatabaseManifest = filepath.Join(SeqDatabaseDir, "manifest.json")
)

var (
	// DefaultConfig is the initiate client config that's embedded with repp
	// and installed on the first run
	//go:embed config.yaml
	DefaultConfig []byte

	// DefaultEnzymes is the JSON file of default enzymes embedded with repp
	//go:embed enzymes.json
	DefaultEnzymes []byte

	// DefaultFeatures is the JSON file of default features embedded with repp
	//go:embed features.json
	DefaultFeatures []byte

	// DefaultPrimer3Config is the FS of Primer3, needed to run primer3_core, etc
	//go:embed primer3_config primer3_config/interpretations
	DefaultPrimer3Config embed.FS
)

// SynthCost contains data of the cost of synthesizing DNA up to a certain
// size. Can be fixed (ie everything beneath that limit is the same amount)
// or not (pay by the bp)
type SynthCost struct {
	// whether it's a fixed or variable cost
	Fixed bool `mapstructure:"fixed"`

	// the cost (either per bp or for the whole stretch)
	Cost float64 `mapstructure:"cost"`
}

// Config is the Root-level settings struct and is a mix
// of settings available in config.yaml and those
// available from the command line
type Config struct {
	// the config file's version
	Version string `mapstructure:"version"`

	// the cost of each Gibson Assembly
	GibsonAssemblyCost float64 `mapstructure:"gibson-assembly-cost"`

	// the cost of time for each Gibson Assembly
	GibsonAssemblyTimeCost float64 `mapstructure:"gibson-assembly-time-cost"`

	// the cost per bp of synthesized DNA as a fragment (as a step function)
	SyntheticFragmentCost map[int]SynthCost `mapstructure:"synthetic-fragment-cost"`

	// the cost per bp of synthesized clonal DNA  (delivered in a plasmid)
	SyntheticPlasmidCost map[int]SynthCost `mapstructure:"synthetic-plasmid-cost"`

	// the maximum number of fragments in the final assembly
	FragmentsMaxCount int `mapstructure:"fragments-max-count"`

	// the minimum homology between this fragment and the net one
	FragmentsMinHomology int `mapstructure:"fragments-min-junction-length"`

	// maximum length of homology between two adjacent fragments in bp
	FragmentsMaxHomology int `mapstructure:"fragments-max-junction-length"`

	// maximum allowable hairpin melting temperature (celcius)
	FragmentsMaxHairpinMelt float64 `mapstructure:"fragments-max-junction-hairpin"`

	// the cost per bp of primer DNA
	PcrBpCost float64 `mapstructure:"pcr-bp-cost"`

	// the cost of each PCR reaction
	PcrRxnCost float64 `mapstructure:"pcr-rxn-cost"`

	// the cost of time for each PCR reaction
	PcrTimeCost float64 `mapstructure:"pcr-time-cost"`

	// PcrMinLength is the minimum size of a fragment (used to filter BLAST results)
	PcrMinLength int `mapstructure:"pcr-min-length"`

	// the maximum primer3 score allowable
	PcrPrimerMaxPairPenalty float64 `mapstructure:"pcr-primer-max-pair-penalty"`

	// the maximum length of a sequence to embed up or downstream of an amplified sequence
	PcrPrimerMaxEmbedLength int `mapstructure:"pcr-primer-max-embed-length"`

	// PcrPrimerMaxOfftargetTm is the maximum tm of an offtarget, above which PCR is abandoned
	PcrPrimerMaxOfftargetTm float64 `mapstructure:"pcr-primer-max-ectopic-tm"`

	// PcrBufferLength is the length of buffer from the ends of a match in which
	// to allow Primer3 to look for a primer
	PcrBufferLength int `mapstructure:"pcr-buffer-length"`

	// minimum length of a synthesized piece of DNA
	SyntheticMinLength int `mapstructure:"synthetic-min-length"`

	// maximum length of a synthesized piece of DNA
	SyntheticMaxLength int `mapstructure:"synthetic-max-length"`
}

// Setup checks that the REPP data directory exists.
// It creates one and writes default config files to it otherwise.
func Setup() {
	// create the REPP directory if it doesn't exist
	_, err := os.Stat(reppDir)
	if os.IsNotExist(err) {
		if err = os.Mkdir(reppDir, 0755); err != nil {
			log.Fatal(err)
		}
	} else if err != nil {
		log.Fatal(err)
	}

	// create the sequence database directory if it doesn't exist
	_, err = os.Stat(SeqDatabaseDir)
	if os.IsNotExist(err) {
		if err = os.Mkdir(SeqDatabaseDir, 0755); err != nil {
			log.Fatal(err)
		}
	} else if err != nil {
		log.Fatal(err)
	}

	// copy the default config file if it doesn't exist
	_, err = os.Stat(configPath)
	if os.IsNotExist(err) {
		if err = os.WriteFile(configPath, DefaultConfig, 0644); err != nil {
			log.Fatal(err)
		}
	} else if err != nil {
		log.Fatal(err)
	}

	// create the features DB if it doesn't exist
	_, err = os.Stat(FeatureDB)
	if os.IsNotExist(err) {
		if err = os.WriteFile(FeatureDB, DefaultFeatures, 0644); err != nil {
			log.Fatal(err)
		}
	} else if err != nil {
		log.Fatal(err)
	}

	// create the enzymes DB if it doesn't exist
	_, err = os.Stat(EnzymeDB)
	if os.IsNotExist(err) {
		if err = os.WriteFile(EnzymeDB, DefaultEnzymes, 0644); err != nil {
			log.Fatal(err)
		}
	} else if err != nil {
		log.Fatal(err)
	}

	// create the primer3 config directory if it does not exist
	_, err = os.Stat(Primer3Config)
	if os.IsNotExist(err) {
		copyFromEmbeded(DefaultPrimer3Config, "primer3_config", Primer3Config)
	} else if err != nil {
		log.Fatal(err)
	}
}

// copyFrom copies an embedded directory to a local directory recursively
func copyFromEmbeded(fs embed.FS, from, to string) {
	if err := os.Mkdir(to, 0755); err != nil && !os.IsExist(err) {
		log.Fatal(err)
	}

	entries, err := fs.ReadDir(from)
	if err != nil {
		log.Fatal(err)
	}

	for _, entry := range entries {
		if entry.IsDir() {
			copyFromEmbeded(fs, path.Join(from, entry.Name()), path.Join(to, entry.Name()))
			continue
		}

		contents, err := fs.ReadFile(path.Join(from, entry.Name()))
		if err != nil {
			log.Fatal(err)
		}
		if err = os.WriteFile(path.Join(to, entry.Name()), contents, 0644); err != nil {
			log.Fatal(err)
		}
	}
}

// New returns a new Config struct populated by settings from
// config.yaml, in the repo, or some other settings file the user
// points to with the "--config" command
//
// TODO: check for and error out on nonsense config values
// TODO: add back the config file path setting
func New() *Config {
	// read in the default settings first
	viper.SetConfigType("yaml")
	viper.SetConfigFile(configPath)
	if err := viper.ReadInConfig(); err != nil {
		log.Fatal(err)
	}

	if userConfig := viper.GetString("config"); userConfig != "" {
		viper.SetConfigFile(userConfig)               // user has specified a new path for a settings file
		if err := viper.MergeInConfig(); err != nil { // read in user defined settings file
			log.Fatal(err)
		}

		file, _ := os.Open(userConfig)
		userData := make(map[string]interface{})
		if err := yaml.NewDecoder(file).Decode(userData); err != nil {
			log.Fatal(err)
		}

		userConfig := &Config{}
		if err := mapstructure.Decode(userData, userConfig); err != nil {
			log.Fatal(err)
		}
	}

	config := &Config{}
	if err := viper.Unmarshal(&config); err != nil {
		log.Fatalf("failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}
	return config
}

// SynthFragmentCost returns the cost of synthesizing a linear stretch of DNA
func (c *Config) SynthFragmentCost(fragLength int) float64 {
	// by default, we try to synthesize the whole thing in one piece
	// we may optionally need to split it into multiple
	fragCount := math.Ceil(float64(fragLength) / float64(c.SyntheticMaxLength))
	fragLength = int(math.Floor(float64(fragLength) / float64(fragCount)))

	cost := synthCost(fragLength, c.SyntheticFragmentCost)
	if cost.Fixed {
		return fragCount * cost.Cost
	}

	return fragCount * float64(fragLength) * cost.Cost
}

// SynthPlasmidCost returns the cost of synthesizing the insert and having it delivered in a plasmid
func (c *Config) SynthPlasmidCost(insertLength int) float64 {
	cost := synthCost(insertLength, c.SyntheticPlasmidCost)
	if cost.Fixed {
		return cost.Cost
	}

	return float64(insertLength) * cost.Cost
}

// synthCost returns the cost of synthesizing a piece of DNA
func synthCost(seqLength int, costs map[int]SynthCost) SynthCost {
	// find the smallest synth length greater than fragLength
	// Ex: a synthesis provider may say it's 32 cents up to 500bp and
	// 60 cents up to 2000bp. So, for a 750bp sequence, we want to use
	// the 2000bp price
	// TODO: add error here for if there's no cost for seqLength (too large)
	costLengthKeys := make([]int, len(costs))
	for key := range costs {
		costLengthKeys = append(costLengthKeys, key)
	}
	sort.Ints(costLengthKeys)

	synthCostKey := 0
	for _, keyLength := range costLengthKeys {
		if keyLength >= seqLength {
			synthCostKey = keyLength
			break
		}
	}

	// we're not able to make a fragment/gene this large
	// return an extremely large number
	if synthCostKey == 0 {
		return SynthCost{
			Fixed: true,
			Cost:  math.MaxInt32,
		}
	}

	return costs[synthCostKey]
}
