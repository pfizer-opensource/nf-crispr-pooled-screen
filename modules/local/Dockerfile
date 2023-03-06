From python:3.10-buster

Run apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# pip install from requirements.txt
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

