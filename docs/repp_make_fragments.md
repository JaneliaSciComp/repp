---
layout: default
title: fragments
parent: make
grand_parent: repp
nav_order: 2
---
## repp make fragments

Build a plasmid from its constituent fragments

### Synopsis

Prepare a list of fragments for assembly via Gibson Assembly. Fragments are
checked for existing homology with their neighbors and are prepared for
assembly with PCR.

```
repp make fragments [flags]
```

### Options

```
  -b, --backbone string             backbone to insert the fragments into. Can either be an entry 
                                    in one of the dbs or a file on the local filesystem.
  -d, --dbs string                  comma separated list of sequence databases by name
  -e, --enzymes string              comma separated list of enzymes to linearize the backbone with.
                                    The backbone must be specified. 'repp ls enzymes' prints a list of
                                    recognized enzymes.
  -h, --help                        help for fragments
  -i, --in string                   input file name (FASTA or Genbank)
  -o, --out string                  output file name
      --synthetic-frag-factor int   Penalty for synthetic fragments (default 1)
```

### Options inherited from parent commands

```
  -c, --config string           User defined config file that may override all or some default settings
      --primer3-config string   primer3 config folder to be used instead of the default
      --repp-data-dir string    Default REPP data directory
  -v, --verbose                 write DEBUG logs
```

### SEE ALSO

* [repp make](repp_make)	 - Make a plasmid from its expected sequence, features, or fragments

###### Auto generated by spf13/cobra on 21-Sep-2023
