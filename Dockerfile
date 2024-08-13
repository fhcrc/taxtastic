# docker build --build-arg VERSION=something --tag taxtastic:latest .
FROM python:3.11-slim-bullseye

ARG VERSION=
ENV PIP_ROOT_USER_ACTION=ignore TAXTASTIC_VERSION=${VERSION#v}

RUN apt-get -y update && apt-get upgrade -y && apt-get install -y unzip wget

WORKDIR /opt/build

COPY dev/install_pplacer.sh ./
RUN /opt/build/install_pplacer.sh /usr/local

COPY setup.py MANIFEST.in README.rst ./
COPY taxtastic/ ./taxtastic/
RUN pip3 install --upgrade pip && pip3 install .

WORKDIR /opt/run
RUN mkdir -p /app /fh /mnt /run/shm

CMD ["taxit", "-h"]
