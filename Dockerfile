FROM alpine:latest AS build

# Compile NetCDF4
RUN mkdir /tmp/setup
RUN apk --no-cache add g++ gcc m4 make zlib-dev libc-dev curl-dev linux-headers
RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/testing hdf5-dev
WORKDIR /tmp/setup
RUN wget -q https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-4.6.1.tar.gz
RUN tar -xzf netcdf-4.6.1.tar.gz
WORKDIR /tmp/setup/netcdf-4.6.1
RUN ./configure --prefix=/usr/local && make && make install
WORKDIR /
RUN rm -rf /tmp/setup

# Compile Proj4
RUN mkdir /tmp/setup
RUN apk --no-cache add git cmake python3 sqlite sqlite-dev
WORKDIR /tmp/setup
RUN git clone https://github.com/OSGeo/proj.4.git proj
WORKDIR proj
RUN mkdir build
WORKDIR build
RUN cmake /tmp/setup/proj
RUN make && make install
WORKDIR /
RUN rm -rf /tmp/setup

# Compile GDAL
RUN mkdir /tmp/setup
RUN apk --no-cache add autoconf
WORKDIR /tmp/setup
RUN git clone https://github.com/OSGeo/gdal.git
WORKDIR gdal/gdal
RUN ./autogen.sh
RUN ./configure --prefix=/usr/local && make && make install
WORKDIR /
RUN rm -rf /tmp/setup

# Add pip dependencies
RUN apk add --no-cache python3-dev openblas-dev freetype-dev
RUN pip3 install --upgrade pip
RUN pip3 install matplotlib
RUN pip3 install netCDF4
RUN pip3 install scipy
RUN pip3 install haversine
RUN pip3 install wget
RUN pip3 install utm
RUN pip3 install requests
RUN pip3 install tqdm
RUN pip3 install eventlet
RUN pip3 install bs4
RUN pip3 install seaborn
RUN pip3 install pandas
RUN pip3 install pyproj

# Install AdcircPy from pip
RUN pip3 install AdcircPy
RUN python3 - << import AdcircPy
RUN pip freeze > requirements.txt
RUN cat requirements.txt
