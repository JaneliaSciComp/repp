DOCKER_IMAGE=jjtimmons/repp
VERSION=v1.0.0

.PHONY: build
build: fmt lint
	go mod tidy
	go mod vendor
	go build -o ./bin/repp ./cmd/repp

.PHONY: install
install:
	go install ./cmd/repp

.PHONY: image
image:
	docker build \
		-t ${DOCKER_IMAGE}:$(VERSION) \
		-t ${DOCKER_IMAGE}:latest .

image/push: image
	docker push ${DOCKER_IMAGE}:$(VERSION)
	docker push ${DOCKER_IMAGE}:latest

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
	gofmt -l ./cmd ./internal

lint:
	golangci-lint run
