## Build Base ##########################################################

FROM alpine:3.20 AS build-base

RUN apk --no-cache add \
      build-base \
      ca-certificates \
      wget

## Build UDUNITS #######################################################

FROM build-base AS build-udunits

ARG UDUNITS_VERSION=2.2.28

RUN apk --no-cache add \
      expat-dev

RUN mkdir -p /tmp/udunits \
 && cd /tmp/udunits \
 && wget --no-check-certificate \
      https://downloads.unidata.ucar.edu/udunits/${UDUNITS_VERSION}/udunits-${UDUNITS_VERSION}.tar.gz \
 && tar xf udunits-${UDUNITS_VERSION}.tar.gz \
 && cd udunits-${UDUNITS_VERSION} \
 && ./configure --prefix=/tmp/udunits/build \
 && make \
 && make install

## Build R Packages ####################################################

FROM build-base AS build-rpkgs

RUN apk --no-cache add \
      R \
      R-dev \
      R-doc

COPY packages.R /tmp/rpkgs/

RUN mkdir -p /tmp/rpkgs/build \
 && cd /tmp/rpkgs \
 && Rscript packages.R

## Final Image #########################################################

FROM alpine:3.20

RUN apk --no-cache add \
      R \
      bash \
      coreutils \
      libexpat \
      font-inconsolata \
      font-dejavu \
      font-noto \
      font-noto-cjk \
      font-awesome \
      font-noto-extra \
      freetype \
      pandoc

ENV OSFONTDIR=/usr/share/fonts
RUN fc-cache --really-force --verbose

COPY --from=build-udunits /tmp/udunits/build /usr/local
COPY --from=build-rpkgs   /tmp/rpkgs/build   /usr/lib/R/library
