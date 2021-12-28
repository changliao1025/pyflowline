FROM gitpod/workspace-full

RUN sudo apt-get update  &&\
    sudo apt-get install gdal-bin  &&\
    sudo apt-get install -y     libgdal-dev  &&\
    export CPLUS_INCLUDE_PATH=/usr/include/gdal  &&\
    export C_INCLUDE_PATH=/usr/include/gdal