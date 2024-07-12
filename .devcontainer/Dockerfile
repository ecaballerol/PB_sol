# For more information, please refer to https://aka.ms/vscode-docker-python
FROM continuumio/miniconda3

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

WORKDIR /app/opt

# Installing python dependencies pip requirements
RUN conda create --name csi python=3.10 numpy

# #it works
SHELL ["conda", "run", "-n", "csi", "/bin/bash", "-c"]

RUN apt-get update && apt-get install -y \
    git \
    build-essential pkg-config libhdf5-serial-dev\
    && apt-get clean

# Clone the okada4py repository
RUN git clone https://github.com/jolivetr/okada4py.git

# Change directory to the cloned repository
WORKDIR /app/opt/okada4py
# Install okada4py within the Conda environment
#Running well
RUN pip install .
RUN mv /app/opt/okada4py/build/lib*/*.so .

WORKDIR /app/opt
# Clone the okada4py repository
RUN git clone https://github.com/jolivetr/csi.git

ENV PYTHONPATH=/app/opt

#other dependencies of python
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

#Copy the files to make the inversion
COPY . /app



# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
# CMD ["python", "StaSLip.py"]

# Set the default command to run a Python shell or your script
CMD ["/bin/bash", "-c", "conda activate csi && exec bash"]
