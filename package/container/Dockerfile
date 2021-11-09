FROM debian:bullseye

WORKDIR /build
COPY . .

RUN apt-get update
RUN apt-get install -y \
  make \
  g++ \
  zlib1g-dev

RUN make clean
RUN make
RUN cp -p ./pstools /usr/local/bin
RUN /usr/local/bin/pstools version

CMD /usr/local/bin/pstools
