<img src="https://user-images.githubusercontent.com/13923102/72353248-b90d7680-36b1-11ea-8714-3249a887b156.png" width="650" margin="0 auto 10px auto" />

`repp` is a tool for DNA assembly. It takes a target plasmid and finds the least expensive combination of fragments from user and public repositories to create it via Gibson Assembly.

Biologists profit when they can re-use DNA during plasmid design: it enables cheaper designs and faster builds. But parsing through all re-usable DNA is completely infeasible. For example, there are over 75,000 plasmids in Addgene -- the likelihood of knowing the best combination and ordering of sub-sequences from Addgene for a given plasmid design is low.

`repp` does such plasmid design. It turns specifications into assembly plans that use the least expensive combination of existing (PCR) and newly synthesized DNA fragments. Target plasmids are specifiable using their target sequence, features, or sub-fragments.

## Publication

We published a paper about `repp` in PLOS One: [Timmons, J.J. & Densmore D. Repository-based plasmid design. PLOS One.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0223935) We used it to build thousands of plasmids from iGEM and Addgene and showed that it reduced the cost of plasmid design as compared to synthesis.

## Examples

See [/examples](/examples) to see input/output from `repp`.

<br>

![https://user-images.githubusercontent.com/13923102/72355113-d55ee280-36b4-11ea-8663-f5759cd7597b.png](https://user-images.githubusercontent.com/13923102/72355113-d55ee280-36b4-11ea-8663-f5759cd7597b.png)

## Documentation

See [the docs](https://lattice-automation.github.io/repp/) or `--help` for any command.

## Installation

### From Docker

Run `repp` via Docker:

```sh
alias repp="docker run -i --rm -v repp:/root/.repp jjtimmons/repp:latest"
repp --help
```

### From Source

Note: `repp` depends on:

- [`Go >= 1.18.0`](https://go.dev/doc/install) for compilation
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for sequence alignment
- [Primer3](https://github.com/primer3-org/primer3) for hairpin detection and off-target primer-binding detection

```sh
git clone https://github.com/Lattice-Automation/repp.git
cd repp
make install
```

## Sequence Databases

`repp` uses sequence databases for plasmid assembly. These are added as FASTA files along with the name and cost per plasmid from that source.

Some existing FASTA files are maintained in our S3 bucket [`repp`](https://s3.console.aws.amazon.com/s3/buckets/repp?region=us-east-1&tab=objects). Below is a snippet for downloading and installing each via the [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html):

```sh
for db in igem addgene dnasu; do
  aws s3 cp "s3://repp/$db.fa.gz" .
  gzip -d "$db.fa.gz"
done

# add sequence DBs with the cost of ordering a plasmid from each source
repp add database --name igem --cost 0.0 < igem.fa
repp add database --name addgene --cost 65.0 < addgene.fa
repp add database --name dnasu --cost 55.0 < dnasu.fa
```

## Plasmid Design

### Sequence

To design a plasmid based on its target sequence save it to a FASTA file. For example:

```
>2ndVal_mScarlet-I
CAACCTTACCAGAGGGCGCCCCAGCTGGCAATTCCGACGTCTAAGAAACCATTATTATCA...
```

Then call `repp make sequence` to design it. The following example uses Addgene and a local BLAST database `parts_library.fa` as fragment sources:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --dbs "addgene,parts_library.fa"
```

### Features

To design a plasmid based on the features it should contain, specify the features by name. By default, these should refer to features that are in `repp`'s feature database (`~/.repp/features.tsv`). Features can also refer to fragments, as in the following example where a plasmid is specified by its constituent list of iGEM parts:

```bash
repp make features "BBa_R0062,BBa_B0034,BBa_C0040,BBa_B0010,BBa_B0012" --backbone pSB1C3 --enzymes "EcoRI,PstI" --dbs igem
```

### Fragments

To design a plasmid from its constituent fragments, save them to a multi-FASTA.

```txt
>GFP
ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGG...
>backbone
TACTAGTAGCGGCCGCTGCAGTCCGGCAAAAAAGGGCAAGGTGTCACCACCCTGCCCTT...
```

And call the file from `repp make fragments`:

```bash
repp make fragments --in "./fragments.fa" --out "plasmid.json"
```

### Configuration

The default settings file used by `repp` is in `~/.repp/config.yaml`. The maximum number of fragments in an assembly, the minimum overlap between adjacent fragments, and cost curves for synthesis are all defined there. Editing this file directly will change the default values used during plasmid designs. For more details, see [configuration](https://lattice-automation.github.io/repp/configuration).

To overwrite some `repp` settings on a per-design basis, create another YAML file:

```yaml
# custom_settings.yaml
fragments-min-junction-length: 25
synthetic-fragment-cost:
  1800: # max length in bp
    fixed: false
    cost: 0.07 # per bp cost
```

And reference it during plasmid design:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --dbs addgene --settings "./custom_settings.yaml"
```

### Backbones and Enzymes

The plasmid sequence in the input file is designed as a circular plasmid by default. In other words, `repp` assumes that the sequence includes an insert sequence as well as a backbone. To use the sequence in the input file as an insert sequence but another fragment as a backbone, use the `--backbone` and `--enzymes` command in combination. This will lookup `--backbone` in the fragment databases and digest it with the enzyme selected through the `--enzymes` flag. The linearized backbone will be concatenated to the insert sequence. For example, to insert a `GFP_CDS` sequence into iGEM's `pSB1A3` backbone after linearizing it with `PstI` and `EcoRI`:

```bash
repp make sequence --in "./GFP_CDS.fa" --dbs addgene,igem --backbone pSB1A3 --enzymes "PstI,EcoRI"
```

The largest linearized fragment post-digestion with all enzymes is used as the backbone in the Gibson Assembly.

### Output

`repp` saves plasmid designs to JSON files at the path specified through the `--out` flag. Below is an abbreviated example of plasmid design output:

```json
{
  "target": "2ndVal_mScarlet-I",
  "seq": "CAACCTTACCAGAGGGCGCCCCAG...",
  "time": "2019/06/24 11:51:39",
  "solutions": [
    {
      "count": 2,
      "cost": 236.65,
      "fragments": [
        {
          "type": "pcr",
          "cost": 94.67,
          "url": "https://www.addgene.org/103998/",
          "seq": "ACAAATAAATGTCCAGACCTGCAG...",
          "pcrSeq": "ACAAATAAATGTCCAGACCTGCAG...",
          "primers": [
            {
              "seq": "ACAAATAAATGTCCAGACCTGCA",
              "strand": true,
              "penalty": 4.406456,
              "pairPenalty": 18.046225,
              "tm": 58.594,
              "gc": 39.13
            },
            {
              "seq": "CATATGTATATCTCCTTCTTAAATCT",
              "strand": false,
              "penalty": 13.639769,
              "pairPenalty": 18.046225,
              "tm": 52.36,
              "gc": 26.923
            }
          ]
        },
        {
          "id": "103998-103998-synthesis-1",
          "type": "synthetic",
          "cost": 129,
          "seq": "AGGAGATATACATATGGTGAGTAA..."
        }
      ]
    }
  ]
}
```
