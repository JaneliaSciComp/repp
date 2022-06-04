.PHONY: build
build: fmt
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
	cd docs && make

docs/serve: docs
	cd docs && make serve

fmt:
	gofmt -l ./internal

lint:
	golangci-lint run
