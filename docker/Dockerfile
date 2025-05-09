FROM ubuntu:25.04

# Avoid interactive prompts (for tzdata, etc.)
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
  python3.13 python3-pip openjdk-17-jre-headless procps curl \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Optionally: Symlink python3 → python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install uv as you do before
RUN curl -LsSf https://astral.sh/uv/install.sh | sh


# Set up working directory
WORKDIR /app

# Copy project files
COPY . /app

# Make sure bcftools is executable
RUN chmod +x /app/tools/bcftools

# Install VCFstash with uv
RUN uv venv /app/.venv && \
    . /app/.venv/bin/activate && \
    uv pip install hatchling editables && \
    uv pip install -e ".[parquet]"

# Set environment variable for the application
ENV VCFSTASH_ROOT=/app
ENV PATH="/app/.venv/bin:${PATH}"

# Set the entry point to the vcfstash CLI
ENTRYPOINT ["vcfstash"]