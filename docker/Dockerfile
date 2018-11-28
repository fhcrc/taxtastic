FROM python:3.6.7-stretch

COPY taxtastic.tar.gz /src/taxtastic.tar.gz
COPY install_pplacer.sh /usr/local/bin/install_pplacer.sh

RUN apt-get update && \
    apt-get install -y unzip --no-install-recommends && \
    /usr/local/bin/install_pplacer.sh /usr/local 1.1.alpha19 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get purge -y --auto-remove unzip && \
    pip install /src/taxtastic.tar.gz

RUN mkdir -p /app /fh /mnt /run/shm

CMD ["taxit", "-h"]
