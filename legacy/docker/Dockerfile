# Use Ubuntu 18.04 as the base image
FROM ubuntu:18.04
WORKDIR /IMBAtracker
# Set noninteractive mode and timezone
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
# Update package lists and install necessary dependencies
RUN apt-get update && \
    apt-get install -y --fix-missing build-essential zlib1g-dev \
    libncurses5-dev libgdbm-dev libnss3-dev \
    libssl-dev libreadline-dev libffi-dev curl software-properties-common \pyqt5-dev-tools \
    qttools5-dev-tools 

RUN apt-get install -y  --fix-missing\
    liblpsolve55-dev \
    libalglib-dev \
    libboost-all-dev \
    git unzip 

RUN apt-get install -y  --fix-missing libopencv-dev

RUN apt-get install -y cmake 

WORKDIR /app

# Copy your project files into the container
COPY . /app
# Change directory to IMBtracker_cpp
WORKDIR /app/IMBtracker_cpp

RUN unzip cvblob.zip


WORKDIR /app/IMBtracker_cpp/cvblob
RUN cmake .
RUN make

WORKDIR /app/IMBtracker_cpp
RUN cmake .
RUN make

# CMD specifies the command to run your Python script when the container starts
CMD ["./lrvTrack"]

# ENV QT_QPA_PLATFORM=offscreen
# # Update package lists and install wget and tar
# RUN apt-get update && \
#     apt-get install -yqq --no-install-recommends wget tar

# # Download and extract Python 3.9 source code
# RUN wget https://www.python.org/ftp/python/3.9.0/Python-3.9.0.tar.xz && \
#     tar -xf Python-3.9.0.tar.xz && \
#     cd Python-3.9.0 && \
#     ./configure && \
#     make && \
#     make install

# # Install pip and setuptools
# RUN apt-get install -y python3-pip python3-setuptools python3-pyqt5

# # Install additional Python dependencies
# RUN pip3 install --upgrade pip scikit-build opencv-contrib-python-headless 
# RUN pip install PyQt5
# # Set the working directory

