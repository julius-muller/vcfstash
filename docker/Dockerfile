# docker/Dockerfile
FROM python:3.11-slim

# Install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    openjdk-11-jre-headless \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && curl -LsSf https://astral.sh/uv/install.sh | sh

# Set up working directory
COPY . /app
WORKDIR /app

# Install Python dependencies
RUN uv venv && uv pip install -r requirements.txt && uv pip install .

# Make sure bcftools is executable
RUN chmod +x /app/tools/bcftools

# Set environment variable for the application
ENV VCFSTASH_ROOT=/app

# Set the entry point
ENTRYPOINT ["python", "/app/vcfstash.py"]
