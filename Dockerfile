FROM golang:1.18.3

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ncbi-blast+ \
        primer3

WORKDIR $HOME/src
ADD . .
RUN go install ./cmd/repp
ENTRYPOINT ["repp"]
