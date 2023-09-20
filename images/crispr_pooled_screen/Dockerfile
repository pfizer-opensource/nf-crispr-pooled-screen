FROM python:3.10-buster

RUN apt update \
 && apt install --assume-yes \
      build-essential \
 && rm -rf /var/lib/apt/lists/*

# pip install from requirements.txt
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir \
      --requirement /tmp/requirements.txt \
 && rm /tmp/requirements.txt
