FROM gitpod/workspace-full-vnc:latest

# Qt5 graphics libraries for napari
RUN sudo apt-get update && \
    sudo apt-get install -y \
        libpython3.8-dev \
        qtbase5-dev \
        qtchooser \
        qt5-qmake \
        qtbase5-dev-tools && \
    sudo rm -rf /var/lib/apt/lists/*

RUN python -m pip install "napari[all]"

# Install cylinter
COPY --chown=gitpod:gitpod . /app
RUN python -m pip install /app

# Client for synapse.org
RUN python -m pip install synapseclient
