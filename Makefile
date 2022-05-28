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

.PHONY: dist
dist:
	go mod tidy
	go build -o ./bin/repp
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

all: build install

dbs:
	cd assets && sh makeblastdbs.sh

uninstall: clean
	rm $(APP)
	rm -rf $(APP_DATA)

test: all
	go test -timeout 200s ./internal/repp

.PHONY: docs
docs:
	go run ./docs/main.go
	cp README.md ./docs/index.md
	cd docs && bundle exec just-the-docs rake search:init
	find ./docs -name *make* -type f -exec sed -i -e 's/\/Users\/josh/~/g' {} \;
	rm ./docs/*-e

docs-dev: docs 
	cd docs && bundle exec jekyll serve