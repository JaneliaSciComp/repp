NAME=repp

.DEFAULT_GOAL: build

.PHONY: dist
build:
	go mod tidy && \
		go mod vendor && \
		go build -o ./bin/repp ./cmd

install:
	mkdir -p $(APP_DATA)

	cp $(SETTINGS) $(APP_DATA)/config.yaml
	cp -r ./vendor/primer3_config $(APP_DATA) 
	cp -r ./assets/addgene/db/** $(APP_DATA) 
	cp -r ./assets/igem/db/** $(APP_DATA)
	cp -r ./assets/dnasu/db/** $(APP_DATA)
	cp ./assets/snapgene/features.tsv $(APP_DATA)
	cp ./assets/neb/enzymes.tsv $(APP_DATA)

.PHONY: test
test: all
	go test -timeout 200s ./internal/repp

.PHONY: docs
docs:
	go run ./docs/main.go
	cp README.md ./docs/index.md
	cd docs && bundle exec just-the-docs rake search:init
	find ./docs -name *make* -type f -exec sed -i -e 's/\/Users\/josh/~/g' {} \;
	rm ./docs/*-e

serve/docs: docs 
	cd docs && bundle exec jekyll serve