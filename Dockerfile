FROM ubuntu:14.04

MAINTAINER Tomas Tinoco De Rubira <tomast@eeh.ee.ethz.ch>

ENV PFNET /packages/pfnet

RUN apt-get update -y
RUN apt-get install -y build-essential
RUN apt-get install -y python-numpy
RUN apt-get install -y python-scipy
RUN apt-get install -y libgraphviz-dev
RUN apt-get install -y cython
RUN apt-get install -y python-nose
RUN apt-get install -y libmumps-dev
RUN apt-get install -y libmumps-seq-dev
RUN apt-get install -y python-matplotlib
RUN apt-get install -y x11-apps
RUN apt-get install -y gpicview

ADD ./packages/pfnet.tar.gz /packages/
ADD ./packages/optalg.tar.gz /packages/

RUN cd /packages/pfnet; rm ./data/*.raw; make lib NO_RAW_PARSER=1
ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PFNET}/lib;/usr/local/lib
RUN cd /packages/pfnet/python; python setup.py install --no_raw_parser; nosetests -v -s

RUN cd /packages/optalg; python setup.py install

ADD ./ /gridopt/
RUN cd /gridopt; rm ./tests/resources/*.raw; python setup.py install; nosetests -v -s

ENTRYPOINT ["gridopt"]

CMD ["--help"]

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
