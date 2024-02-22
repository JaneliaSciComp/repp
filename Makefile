DOCKER_IMAGE=janeliascicomp/repp
DOCKER_PLATFORM_ARG=--platform linux/x86_64,linux/arm64
VERSION=v1.0.0

.PHONY: all
all: build test

ifeq ($(OS), Windows_NT)
    REPP_EXECUTABLE=repp.exe
    RM=rd /S /Q
else
    REPP_EXECUTABLE=repp
    RM=rm -rf
endif

.PHONY: modinstall
modinstall:
	go mod tidy
	go mod vendor

check-everything-committed:
	git diff-index --quiet HEAD

generate-version: check-everything-committed
	git rev-parse HEAD > internal/cmd/commit.txt

.PHONY: build
build: modinstall fmt lint
	go build -o ./bin/$(REPP_EXECUTABLE) ./cmd/repp

build-release: generate-version build

.PHONY: install
install:
	go install ./cmd/repp

.PHONY: image/multiplatform
image/multiplatform:
	docker buildx build ${DOCKER_PLATFORM_ARG} \
		-t ${DOCKER_IMAGE}:$(VERSION) \
		-t ${DOCKER_IMAGE}:latest \
		--push .

release: image/multiplatform
	gh release create $(VERSION) -t $(VERSION) --generate-notes -d

.PHONY: test
test:
	go test $(VERBOSE:%=-v) -timeout 200s ./internal/...

.PHONY: docs
docs:
	cd docs && rm *.md && make

docs/serve: docs
	cd docs && make serve

fmt:
	gofmt -l ./cmd ./internal

lint:
	golangci-lint run

clean:
	go clean -cache -testcache
	$(RM) bin vendor
