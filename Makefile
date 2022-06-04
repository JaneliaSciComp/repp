LOCAL_BIN=/usr/local/bin
APP=${LOCAL_BIN}/repp
APP_DATA=$${HOME}/.repp
SETTINGS=./config/config.yaml

NAME=repp
VERSION=0.1.0

DIST_WIN_ZIP=${NAME}_windows_${VERSION}.zip
DIST_SRC=${NAME}_src_${VERSION}
DIST_SRC_TAR=${DIST_SRC}.tar.gz

PLATFORM:=$(shell uname)

.PHONY: test dist docs
.DEFAULT_GOAL: build

ifeq ($(PLATFORM),Windows_NT)
	$(error Windows not supported via make)
endif

build:
	go get -d
	env GOOS=linux go build -o ./bin/linux -v
	env GOOS=darwin go build -o ./bin/darwin -v
	env GOOS=windows go build -o ./bin/repp.exe -v

install:
	mkdir -p $(APP_DATA)

	cp $(SETTINGS) $(APP_DATA)/config.yaml
	cp -r ./vendor/primer3_config $(APP_DATA) 
	cp -r ./assets/addgene/db/** $(APP_DATA) 
	cp -r ./assets/igem/db/** $(APP_DATA)
	cp -r ./assets/dnasu/db/** $(APP_DATA)
	cp ./assets/snapgene/features.tsv $(APP_DATA)
	cp ./assets/neb/enzymes.tsv $(APP_DATA)

ifeq ($(PLATFORM),Linux)
	install ./bin/linux $(APP)
	install -C ./vendor/linux/blastn $(LOCAL_BIN)
	install -C ./vendor/linux/ntthal $(LOCAL_BIN)
	install -C ./vendor/linux/primer3_core $(LOCAL_BIN)
	install -C ./vendor/linux/blastdbcmd $(LOCAL_BIN)
endif

ifeq ($(PLATFORM),Darwin)
	install ./bin/darwin $(APP)
	install -C ./vendor/darwin/blastn $(LOCAL_BIN)
	install -C ./vendor/darwin/ntthal $(LOCAL_BIN)
	install -C ./vendor/darwin/primer3_core $(LOCAL_BIN)
	install -C ./vendor/darwin/blastdbcmd $(LOCAL_BIN)
endif

windows: build
	cd scripts && makensis windows_installer.nsi

all: build install

dbs:
	cd assets && sh makeblastdbs.sh

uninstall: clean
	rm $(APP)
	rm -rf $(APP_DATA)

test: all
	go test -timeout 200s ./internal/repp

dist-dir:
	mkdir -p ${DIST_SRC}
	rsync -r --delete\
	 --exclude={'.git','dist','test','scripts','bin/repp_install.exe','bin/repp.exe','vendor/windows','assets/addgene/addgene.json','assets/dnasu/DNASU*','assets/igem/xml*','assets/neb/*/'}\
	 . ${DIST_SRC}
	tar -czf ${DIST_SRC_TAR} ${DIST_SRC}
	rm -rf ${DIST_SRC}

dist: windows dist-dir
	cp ./README.md ./docs/index.md

	zip ${DIST_WIN_ZIP} ./bin/repp_install.exe

	scp ${DIST_SRC_TAR} jjtimmons@frs.sourceforge.net:/home/frs/project/repplasmid/
	scp ${DIST_WIN_ZIP} jjtimmons@frs.sourceforge.net:/home/frs/project/repplasmid/

	rm ${DIST_SRC_TAR}
	rm ${DIST_WIN_ZIP}

docs:
	go run ./docs/main.go
	cp README.md ./docs/index.md
	cd docs && bundle exec just-the-docs rake search:init
	find ./docs -name *make* -type f -exec sed -i -e 's/\/Users\/josh/~/g' {} \;
	rm ./docs/*-e

docs-dev: docs 
	cd docs && bundle exec jekyll serve