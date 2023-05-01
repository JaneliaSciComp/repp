## Install Windows Dependencies

These notes are based on scoop - as a package manager to install all the
required dependencies on Windows

Also we assume that you use powershell to run all the command line commands only because scoop instructions are based on powershell. However if the paths are setup correctly you should be able to use CMD too.

### Install scoop

Detailed instructions to install scoop are available at [scoop.sh](https://scoop.sh/)
```
irm get.scoop.sh | iex
```

### Install build tools

scoop install git
scoop install make
scoop install gcc
scoop install go
scoop install golangci-lint

### Create a working folder

```
mkdir <working_dir>
```

### Install primer3
```
cd <workking_dir>
mkdir primer3_latest
cd primer3_latest
git clone https://github.com/primer3-org/primer3.git .
cd src
make TESTOPTS=--windows
```

### Install blast tools
go to https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ and download and run the windows installer.
Typically NCBI blast tools will be installed in `C:\Program Files\NCBI\blast-2.14.0+`

### Install repp

```
cd <working_dir>
mkdir repp_latest
cd repp_latest
git clone https://github.com/JaneliaSciComp/repp.git .
.\setenv.ps1
make build
```

As a note setenv.ps1 assumes that blast was installed in the default location and primer3 tools was installed directly under the user's home directory: `C:\Users\<username>` in the `primer3-latest` build sub-directory.

### Run repp
If you don't setup primer3 and blast tools in the global path you will have to run setenv.ps1 first time you open a powershell window to run repp.

```
cd <working_dir>
.\setenv.ps1
bin\repp.exe -h
```
