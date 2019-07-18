---
layout: default
title: Installation
nav_order: 1
permalink: /installation
---

# Installation

Download links are available at SourceForge: [https://sourceforge.net/projects/repplasmid/files/](https://sourceforge.net/projects/repplasmid/files/).

The BLAST database files are too large for Github but not SourceForge. The code alone is at Github without the databases in the `/assets` directory.

## MacOS/Linux

```bash
wget -O repp_src_0.1.0.tar.gz 'https://sourceforge.net/projects/repplasmid/files/repp_src_0.1.0.tar.gz/download'
tar xzf repp_src_0.1.0.tar.gz
cd repp_src_0.1.0
make install
```

## Windows

1. Download the most recent `repp_windows.*.zip` from [SourceForge](https://sourceforge.net/projects/repplasmid/files/)
2. Unzip
3. Run `repp_install.exe`

## SEE ALSO

- [repp](repp) - REPP
