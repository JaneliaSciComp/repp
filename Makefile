.PHONY: build
build:
	go mod tidy
	go mod vendor
	go build -o ./bin/repp ./cmd/repp

.PHONY: install
install:
	go install ./cmd/repp

.PHONY: test
test:
	go test -timeout 200s ./internal/...

.PHONY: docs
docs:
	cd docs
	go run ./main.go
	cp ../README.md ./index.md
	bundle exec just-the-docs rake search:init
	find . -name *make* -type f -exec sed -i -e 's/\/Users\/josh/~/g' {} \;
	rm ./*-e

serve/docs: docs 
	bundle exec jekyll serve