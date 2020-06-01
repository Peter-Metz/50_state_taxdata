FROM rocker/rstudio:4.0.0
FROM rocker/tidyverse

MAINTAINER Peter Metz <pmetzdc@gmail.com>

WORKDIR /ipopt_df

COPY ./r /ipopt_df/r

# IPoptr envrionment directory
ENV IPOPTR_DIR=/ipopt_df/CoinIpopt/build/Ipopt/contrib

# Give the user root access
RUN echo "$USER ALL = NOPASSWD: ALL" >> /etc/sudoers && \

    apt-get update && apt-get install -y \
        gcc g++ gfortran git patch wget pkg-config liblapack-dev libmetis-dev && \

##################
# Download ipoptr source
# http://www.coin-or.org/Ipopt/documentation/node10.html
##################
    cd /ipopt_df && \
    wget http://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.13.tgz && \
    gunzip Ipopt-3.12.13.tgz && \
    tar -xvf Ipopt-3.12.13.tar && \
    rm -rf Ipopt-3.12.13.tar && \
    mv Ipopt-3.12.13 CoinIpopt

COPY ./coinhsl-2019.05.21 /ipopt_df/CoinIpopt/ThirdParty/HSL

RUN cd /ipopt_df/CoinIpopt && \

# Downloading third party solvers
    cd ThirdParty/Blas && \
        ./get.Blas && \
    cd ../Lapack && \
        ./get.Lapack && \
    cd ../ASL && \
        ./get.ASL && \
    cd ../Mumps && \
        ./get.Mumps && \
    cd ../Metis && \
        ./get.Metis && \
    cd ../HSL && \
        ./configure && \

##################
# Compile ipoptr
##################
    cd ../../ && \
    mkdir build && \
    cd build && \
    ../configure -with-pic CXXFLAGS="-fopenmp" FCFLAGS="-fopenmp" CFLAGS="-fopenmp" ADD_FFLAGS=-fPIC ADD_CFLAGS=-fPIC ADD_CXXFLAGS=-fPIC && \
    make -j3 && \
    make test && \
    make install && \

##################
# Pre-install ipoptr
##################
    echo "install.packages('/ipopt_df/CoinIpopt/build/Ipopt/contrib/RInterface', repos=NULL, type='source')" \
    >> installIpoptrPackage.R && \
    r installIpoptrPackage.R && \
    rm installIpoptrPackage.R

# Default CMD from 
# https://github.com/rocker-org/rocker/blob/master/rstudio/Dockerfile
# CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]