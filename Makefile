.PHONY: build
build:
	go mod tidy && \
		go mod vendor && \
		go build -o ./bin/repp ./cmd

.PHONY: install
install:
	go install ./cmd/repp

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