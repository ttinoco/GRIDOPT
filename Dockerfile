FROM jupyterhub/jupyterhub

MAINTAINER Tomas Tinoco De Rubira <tomast@eeh.ee.ethz.ch>

ENV PFNET /packages/pfnet

RUN apt-get update -y; apt-get install -y build-essential python-numpy python-scipy libgraphviz-dev cython python-nose libmumps-dev libmumps-seq-dev python-matplotlib git

RUN apt-get update -y; apt-get install -y python3-matplotlib

RUN pip install cython numpy notebook scipy

# Make PFNET
RUN mkdir /packages
WORKDIR /packages
RUN git clone https://github.com/ttinoco/PFNET.git
WORKDIR /packages/PFNET/
#RUN rm ./data/*.raw
ENV PFNET /packages/PFNET
ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PFNET}/lib:/usr/local/lib
RUN make clean
RUN make NO_RAW_PARSER=1
RUN ln -s /packages/PFNET/lib/libpfnet.so  /usr/lib/libpfnet.so

# Install PFNET python bindings
WORKDIR /packages/PFNET/python 
RUN python2.7 setup.py install --no_raw_parser
RUN python3 setup.py install --no_raw_parser

#RUN nosetests -v -s

# Install OPTALG
WORKDIR /packages
RUN git clone -b python-3 https://github.com/ttinoco/OPTALG.git
WORKDIR /packages/OPTALG
RUN python2.7 setup.py install
RUN python3 setup.py install

# Install GRIDOPT
WORKDIR /packages
RUN git clone -b python-3 https://github.com/ttinoco/GRIDOPT.git
WORKDIR /packages/GRIDOPT
#RUN rm ./tests/resources/*.raw
RUN python2.7 setup.py install
RUN python3 setup.py install

RUN nosetests -v -s

RUN mkdir /notebooks/
ADD PFNET.ipynb /notebooks/PFNET.ipynb

ADD jupyterhub_config.py /srv/jupyterhub/jupyterhub_config.py
# Dockerized GRIDOPT
#-------------------
#
# Build
# docker build -t gridopt_tag .
#
# Interactive shell
# docker run -i -t --entrypoint=/bin/bash gridopt_tag
#
# Executable
# docker run -t gridopt_tag filename
#
# Save image
# docker save -o gridopt_tag.tar gridopt_tag
#
# Load image
# docker load -i gridopt_tag.tar
#
# Delete containers
# docker rm -f $(docker ps -a -q)
#
# Delete images
# docker rmi -f $(docker images -q)
#
# NOTE for windows: start path with // (e.g. //bin/bash)
