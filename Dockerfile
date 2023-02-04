FROM python:3.7.5-slim

# install python dependencies
RUN python -m pip install \
        pandas==1.1.4 \
        numpy==1.19.4 \
        pysam==0.16.0.1 \
        scipy==1.5.4 \
        argparse

# clone the latest version of intronomer from GitHub
RUN apt-get -y update
RUN apt-get -y install git
RUN git clone https://github.com/pdxgx/intronomer intronomer

# set the command
ENTRYPOINT ["python", "/intronomer/intronomer/intronomer.py"]