FROM golang:1.20

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ncbi-blast+ \
        primer3

VOLUME /data/repp
ENV REPP_DATA_DIR=/data/repp

WORKDIR $HOME/src

ADD . .
RUN go install ./cmd/repp
CMD ["repp"]
