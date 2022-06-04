<img src="https://user-images.githubusercontent.com/13923102/72353248-b90d7680-36b1-11ea-8714-3249a887b156.png" width="650" margin="0 auto 10px auto" />

`REPP` is a tool for DNA assembly. It takes a target plasmid and finds the least expensive combination of fragments from user and public repositories to create it via Gibson Assembly.

Biologists profit when they can re-use DNA during plasmid design: it enables cheaper designs and faster builds. But parsing through all re-usable DNA is completely infeasible. For example, there are over 75,000 plasmids in Addgene -- the likelihood of knowing the best combination and ordering of sub-sequences from Addgene for a given plasmid design is low.

`REPP` enables such plasmid design. It turns plasmid specifications into designs using the least expensive design with both existing DNA fragments (PCR) and newly synthesized DNA fragments. Plasmids are specifiable using their target sequence, features, or sub-fragments.

## Publication

We published a paper about REPP in PLOS One: [Timmons, J.J. & Densmore D. Repository-based plasmid design. PLOS One.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0223935) We used it to build thousands of plasmids from iGEM and Addgene and showed that it reduced the cost of plasmid design as compared to synthesis.

## Installation

Download links are available at SourceForge: [https://sourceforge.net/projects/repplasmid/files/](https://sourceforge.net/projects/repplasmid/files/)

### MacOS/Linux

```bash
wget -O repp_src_0.1.0.tar.gz 'https://sourceforge.net/projects/repplasmid/files/repp_src_0.1.0.tar.gz/download'
tar xzf repp_src_0.1.0.tar.gz
cd repp_src_0.1.0
make install
```

### Windows

1. Download the most recent `repp_windows.*.zip` from [SourceForge](https://sourceforge.net/projects/repplasmid/files/)
2. Unzip
3. Run `repp_install.exe`

## Documentation

See [the docs](https://lattice-automation.github.io/repp/) or use `--help` on any command.

## Examples

See [/examples](/examples) to see input/output from REPP.

<br>

![https://user-images.githubusercontent.com/13923102/72355113-d55ee280-36b4-11ea-8663-f5759cd7597b.png](https://user-images.githubusercontent.com/13923102/72355113-d55ee280-36b4-11ea-8663-f5759cd7597b.png)

## Plasmid Design

### Sequence

To design a plasmid based on its expected sequence save it to a FASTA or Genbank file. For example:

```
>2ndVal_mScarlet-I
CAACCTTACCAGAGGGCGCCCCAGCTGGCAATTCCGACGTCTAAGAAACCATTATTATCA...
```

Then call `repp make sequence` to design it. The following example uses Addgene and a local BLAST database `parts_library.fa` as fragment sources:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --dbs "parts_library.fa"
```

### Features

To design a plasmid based on the features it should contain, specify the features by name. By default, these should refer to features that are in REPP's feature database (`~/.repp/features.tsv`). Features can also refer to fragments, as in the following example where a plasmid is specified by its constituent list of iGEM parts:

```bash
repp make features "BBa_R0062,BBa_B0034,BBa_C0040,BBa_B0010,BBa_B0012" --backbone pSB1C3 --enzymes "EcoRI,PstI" --igem
```

### Fragments

To design a plasmid from its constiuent fragments, save them to a multi-FASTA.

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

### Databases

`REPP` includes three embedded databases from large public repositories: [Addgene](https://www.addgene.org/), [iGEM](http://parts.igem.org/Main_Page), and [DNASU](https://dnasu.org/DNASU/Home.do). Each embedded database and its file path after installation are as follows:

- Addgene, `--addgene`, `~/.repp/addgene`
- DNASU, `--dnasu`, `~/.repp/dnasu`
- iGEM, `--igem`, `~/.repp/igem`

Users can also use their or their lab's fragment databases through the `--dbs` as a list of comma-separated fragment [BLAST databases](https://www.ncbi.nlm.nih.gov/books/NBK279688/). An example of a plasmid design using Addgene, DNASU, and multiple user-defined BLAST repositories is below:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --dnasu --dbs "proteins.fa,backbones.fa"
```

### Configuration

The default settings file used by `REPP` is in `~/.repp/config.yaml`. The maximum number of fragments in an assembly, the minimum overlap between adjacent fragments, and cost curves for synthesis are all defined there. Editing this file directly will change the default values used during plasmid designs. For more details, see [configuration](https://Lattice-Automation.github.io/repp/configuration).

To overwrite some `REPP` settings on a per-design basis, create another YAML file:

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
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --settings "./custom_settings.yaml"
```

### Backbones and Enzymes

The plasmid sequence in the input file is designed as a circular plasmid by default. In other words, REPP assumes that the sequence includes an insert sequence as well as a backbone. To use the sequence in the input file as an insert sequence but another fragment as a backbone, use the `--backbone` and `--enzymes` command in combination. This will lookup `--backbone` in the fragment databases and digest it with the enzyme selected through the `--enzymes` flag. The linearized backbone will be concatenated to the insert sequence. For example, to insert a `GFP_CDS` sequence into iGEM's `pSB1A3` backbone after linearizing it with `PstI` and `EcoRI`:

```bash
repp make sequence --in "./GFP_CDS.fa" --addgene --igem --backbone pSB1A3 --enzymes "PstI,EcoRI"
```

The largest linearized fragment post-digestion with all enzymes is used as the backbone in the Gibson Assembly.

### Output

`REPP` saves plasmid designs to JSON files at the path specified through the `--out` flag. Below is an abbreviated example of plasmid design output:

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
