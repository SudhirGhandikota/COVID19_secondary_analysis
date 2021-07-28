FROM ubuntu:latest
FROM r-base:4.0.3
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools git libcurl4-openssl-dev libssl-dev libxml2-dev

# Copy requirements.R to container current directory
COPY ./requirements.R /tmp/

# Installing the required R packages
RUN Rscript /tmp/requirements.R

COPY ./requirements.txt /tmp/

RUN pip3 install require --no-cache-dir -r /tmp/requirements.txt

#RUN \
   #echo 'alias python="/usr/bin/python3"' >> /root/.bashrc && \
   #echo 'alias pip="/usr/bin/pip3"' >> /root/.bashrc && \
   #source /root/.bashrc
RUN cd "$(dirname $(which python3))" \
					&& ln -s python3 python 
