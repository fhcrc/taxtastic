FROM python:2.7.13-stretch

COPY taxtastic.tar.gz /src/taxtastic.tar.gz

RUN pip install /src/taxtastic.tar.gz

CMD ["taxit", "-h"]
