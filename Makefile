VERSION=v1.0.0

.PHONY: build
build: fmt lint docs
	go mod tidy
	go mod vendor
	go build -o ./bin/repp ./cmd/repp

.PHONY: install
install:
	go install ./cmd/repp

.PHONY: image
image:
	docker build -t jjtimmons/repp:$(VERSION) -t jjtimmons/repp:latest .

image/push: image
	docker push jjtimmons/repp:$(VERSION)
	docker push jjtimmons/repp:latest

release: image/push
	gh release create $(VERSION) -t $(VERSION) --generate-notes -d

.PHONY: test
test:
	go test -timeout 200s ./internal/...

.PHONY: docs
docs:
	cd docs && rm *.md && make

docs/serve: docs
	cd docs && make serve

fmt:
	gofmt -l ./internal

lint:
	golangci-lint run
