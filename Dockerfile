FROM gitpod/workspace-full-vnc:latest

# Qt5 graphics libraries for napari
RUN sudo apt-get update && \
    sudo apt-get install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools && \
    sudo rm -rf /var/lib/apt/lists/*

# Install cylinter
COPY . /app
RUN sudo pip install /app
