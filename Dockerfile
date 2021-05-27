# # taxtastic
#

FROM quay.io/biocontainers/biopython:1.78

RUN mkdir -p /src/taxtastic
ADD setup.py /src/taxtastic/
ADD taxit.py /src/taxtastic/
ADD README.rst /src/taxtastic/
ADD taxtastic /src/taxtastic/taxtastic/
RUN export DEBIAN_FRONTEND=noninteractive
RUN ls -l /src/taxtastic/
RUN pip3 install /src/taxtastic/

ADD pplacer/pplacer-Linux-v1.1.alpha19.zip /src/
WORKDIR /src/
RUN unzip pplacer-Linux-v1.1.alpha19.zip
RUN mv /src/pplacer-Linux-v1.1.alpha19/guppy /usr/local/bin
RUN mv /src/pplacer-Linux-v1.1.alpha19/pplacer /usr/local/bin
RUN mv /src/pplacer-Linux-v1.1.alpha19/rppr /usr/local/bin
WORKDIR /root
RUN rm -rf /src/

CMD ["taxit", "-h"]