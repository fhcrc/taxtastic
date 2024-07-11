FROM python:3.11-slim-bullseye

RUN apt-get -y update && apt-get upgrade -y && apt-get install -y \
unzip \
wget \
# Cirro requirement
psutils

WORKDIR /opt/build
COPY dev/install_pplacer.sh /opt/build/install_pplacer.sh
RUN /opt/build/install_pplacer.sh /usr/local

COPY setup.py MANIFEST.in README.rst requirements.txt /opt/build/
COPY taxtastic /opt/build/taxtastic/
RUN python3 -m pip install --root-user-action=ignore --upgrade pip && \
python3 -m pip install --constraint requirements.txt --root-user-action=ignore .

WORKDIR /opt/run
RUN mkdir -p /app /fh /mnt /run/shm

CMD ["taxit", "-h"]
