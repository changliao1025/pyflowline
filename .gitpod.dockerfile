FROM gitpod/workspace-full

RUN sudo apt-get update  &&\
    sudo apt-get install gdal-bin  &&\
    sudo apt-get install -y     libgdal-dev