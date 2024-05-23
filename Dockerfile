FROM python:3.11-slim-bullseye

RUN apt-get -y update && \
    apt-get install -y unzip wget

WORKDIR /opt/build
COPY dev/install_pplacer.sh /opt/build/install_pplacer.sh
RUN /opt/build/install_pplacer.sh /usr/local

COPY setup.py MANIFEST.in README.rst requirements.txt /opt/build
COPY taxtastic /opt/build/taxtastic
RUN python3 -m pip install --constraint requirements.txt . 

WORKDIR /opt/run
RUN mkdir -p /app /fh /mnt /run/shm

CMD ["taxit", "-h"]
