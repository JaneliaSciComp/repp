---
layout: default
title: fragment
parent: list
grand_parent: repp
nav_order: 1
---
## repp list fragment

List fragments in the databases

### Synopsis

List fragments with a passed name in the specified databases

```
repp list fragment [name] [flags]
```

### Examples

```
  repp list fragment pSB1C3 --dbs igem
```

### Options

```
  -d, --dbs string   comma separated list of sequence databases
  -h, --help         help for fragment
```

### Options inherited from parent commands

```
      --repp-data-dir string   Default REPP data directory
  -v, --verbose                write DEBUG logs
```

### SEE ALSO

* [repp list](repp_list)	 - List any of the things that repp uses to build plasmids

###### Auto generated by spf13/cobra on 28-Mar-2025
