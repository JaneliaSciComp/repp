DOCKER_IMAGE=janeliascicomp/repp
VERSION=v1.0.0

.PHONY: all
all: build test

ifeq ($(OS), Windows_NT)
    REPP_EXECUTABLE=repp.exe
    RM=rd /S /Q
    PLATFORM_ARG=
else
    ARCH := $(shell arch)
	ifeq ($(ARCH), arm64)
		PLATFORM_ARG=--platform linux/x86_64
	else
		PLATFORM_ARG=
	endif
    REPP_EXECUTABLE=repp
    RM=rm -rf
endif

.PHONY: build
build: fmt lint
	go mod tidy
	go mod vendor
	go build -o ./bin/$(REPP_EXECUTABLE) ./cmd/repp

.PHONY: install
install:
	go install ./cmd/repp

.PHONY: image
image:
	docker build ${PLATFORM_ARG} \
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

clean:
	$(RM) bin vendor
