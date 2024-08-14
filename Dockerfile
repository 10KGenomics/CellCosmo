FROM dockerhub.genostack.com:8090/jikai/st-base:py310_r36_zh

RUN conda install star=2.7.10b \
    subread=2.0.1 \
    picard=2.18.17 \
    ucsc-gtftogenepred=447 \
    samtools=1.12 -c bioconda -y

COPY requirements.txt requirements.txt
COPY dist/cell_cosmo-1.0.12.tar.gz cell_cosmo-1.0.12.tar.gz
RUN pip install -r requirements.txt
RUN pip install cell_cosmo-1.0.12.tar.gz
RUN rm cell_cosmo-1.0.12.tar.gz requirements.txt


# docker build -t dockerhub.genostack.com:8090/angs/cellcosmo:1.0.0 .

