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
	"strings"

	"github.com/mitchellh/go-homedir"
	"github.com/mitchellh/mapstructure"
	"github.com/spf13/viper"
	"gopkg.in/yaml.v2"
)

var (
	// reppDir is the root directory where repp settings and database files live
	reppDir string

	// defaultConfigPath is the path to a local/default config file
	defaultConfigPath string

	// defaultPrimer3ConfigDir is the path to a primer3 config directory
	// primer3 is (overly) particular about the trailing slash
	defaultPrimer3ConfigDir string

	// FeatureDB is the path to the features file
	FeatureDB string

	// EnzymeDB is the path to the enzymes file
	EnzymeDB string

	// SeqDatabaseDir is the path to a directory of sequence databases.
	SeqDatabaseDir string

	// SeqDatabaseManifest is the path to the manifest file for the sequence databases.
	SeqDatabaseManifest string
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

	// PcrMinFragLength is the minimum size of a fragment (used to filter BLAST results)
	PcrMinFragLength int `mapstructure:"pcr-min-length"`

	// the maximum primer3 score allowable
	PcrPrimerMaxPairPenalty float64 `mapstructure:"pcr-primer-max-pair-penalty"`

	// the maximum length of a sequence to embed up or downstream of an amplified sequence
	PcrPrimerMaxEmbedLength int `mapstructure:"pcr-primer-max-embed-length"`

	// PcrPrimerMaxOfftargetTm is the maximum tm of an offtarget, above which PCR is abandoned
	PcrPrimerMaxOfftargetTm float64 `mapstructure:"pcr-primer-max-ectopic-tm"`

	// PcrBufferLength is the length of buffer from the ends of a match in which
	// to allow Primer3 to look for a primer
	PcrBufferLength int `mapstructure:"pcr-buffer-length"`

	// Minimum primers length
	PcrPrimerMinLength int `mapstructure:"pcr-min-primer-length"`

	// Maximum primers length
	PcrPrimerMaxLength int `mapstructure:"pcr-max-primer-length"`

	// Optimum primers length
	PcrPrimerOptimumLength int `mapstructure:"pcr-optimum-primer-length"`

	// Min primer annealing temperature (Tm)
	PcrPrimerMinTm float64 `mapstructure:"pcr-primer-min-tm"`

	// Max primer annealing temperature (Tm)
	PcrPrimerMaxTm float64 `mapstructure:"pcr-primer-max-tm"`

	// Max allowed difference in primer annealing temperatures (Tm)
	// If <0 the difference is not checked
	PcrMaxFwdRevPrimerTmDiff float64 `mapstructure:"pcr-max-fwd-rev-primer-tm-diff"`

	// minimum length of a synthesized piece of DNA
	SyntheticMinLength int `mapstructure:"synthetic-min-length"`

	// maximum length of a synthesized piece of DNA
	SyntheticMaxLength int `mapstructure:"synthetic-max-length"`

	// configurable penalty for synthetic fragments
	SyntheticFragmentFactor int `mapstructure:"synthetic-fragment-factor"`

	// user provided path to primer3 config dir
	p3ConfigDir string
}

func initDataPaths(providedReppDir string) (err error) {
	if providedReppDir == "" {
		// if no repp dir was provided
		// try to get it from the environment
		reppDir = os.Getenv("REPP_DATA_DIR")
		if reppDir == "" {
			// use $HOMEDIR/.repp
			var home string
			home, err = homedir.Dir()
			if err != nil {
				return
			}
			reppDir = filepath.Join(home, ".repp")
		}
	} else {
		reppDir = providedReppDir
	}

	defaultConfigPath = filepath.Join(reppDir, "config.yaml")
	defaultPrimer3ConfigDir = filepath.Join(reppDir, "primer3_config") + string(os.PathSeparator)
	FeatureDB = filepath.Join(reppDir, "features.json")
	EnzymeDB = filepath.Join(reppDir, "enzymes.json")
	SeqDatabaseDir = filepath.Join(reppDir, "dbs")
	SeqDatabaseManifest = filepath.Join(SeqDatabaseDir, "manifest.json")

	return
}

// Setup checks that the REPP data directory exists.
// It creates one and writes default config files to it otherwise.
func Setup(providedReppDir string) {

	err := initDataPaths(providedReppDir)
	if err != nil {
		log.Fatal("Error creating repp data paths", err)
	}

	// create the REPP directory if it doesn't exist
	_, err = os.Stat(reppDir)
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

	// the rest of the configuration files are always overwritten for now

	// only copy default config file
	// if it does not exist
	if isConfigFileNeeded(defaultConfigPath, true) {
		if err = os.WriteFile(defaultConfigPath, DefaultConfig, 0644); err != nil {
			log.Fatal(err)
		}
	}

	// features DB
	if isConfigFileNeeded(FeatureDB, true) {
		if err = os.WriteFile(FeatureDB, DefaultFeatures, 0644); err != nil {
			log.Fatal(err)
		}
	}

	// enzymes DB
	if isConfigFileNeeded(EnzymeDB, true) {
		if err = os.WriteFile(EnzymeDB, DefaultEnzymes, 0644); err != nil {
			log.Fatal(err)
		}
	}

	// primer3 config directory
	if isConfigFileNeeded(defaultPrimer3ConfigDir, true) {
		copyFromEmbeded(DefaultPrimer3Config, "primer3_config", defaultPrimer3ConfigDir)
	}
}

func isConfigFileNeeded(configFile string, overwritePreference bool) bool {
	// write the default config file
	if overwritePreference {
		return true
	}
	_, err := os.Stat(configFile)
	if os.IsNotExist(err) {
		return true
	} else if err != nil {
		log.Fatal(err)
	}
	return false
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
	viper.SetConfigFile(defaultConfigPath)
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

// Return the path to the primer3 config directory
func (c *Config) SetPrimer3ConfigDir(p3ConfigDir string) *Config {
	if p3ConfigDir != "" {
		if strings.HasSuffix(p3ConfigDir, "/") {
			c.p3ConfigDir = p3ConfigDir
		} else {
			c.p3ConfigDir = p3ConfigDir + "/"
		}
	}
	return c
}

// Return the path to the primer3 config directory
func (c *Config) GetPrimer3ConfigDir() string {
	if c.p3ConfigDir != "" {
		return c.p3ConfigDir
	} else {
		return defaultPrimer3ConfigDir
	}
}

func (c *Config) SetSyntheticFragmentFactor(value int) *Config {
	c.SyntheticFragmentFactor = value
	return c
}

func (c *Config) GetSyntheticFragmentFactor() int {
	if c.SyntheticFragmentFactor > 0 {
		return c.SyntheticFragmentFactor
	} else {
		return 1
	}
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

func (c *Config) EstimatePCRPrimersLength(defaultValue int) int {
	medPcrPrimerLength := (c.PcrPrimerMinLength + c.PcrPrimerMaxLength) / 2
	if medPcrPrimerLength > 0 {
		return medPcrPrimerLength
	}
	return defaultValue
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
