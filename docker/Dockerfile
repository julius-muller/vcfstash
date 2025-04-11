# docker/Dockerfile
FROM ensemblorg/ensembl-vcf:latest

RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-venv \
    && curl -LsSf https://astral.sh/uv/install.sh | sh

COPY . /app
WORKDIR /app
RUN uv venv && uv pip install .

ENTRYPOINT ["vcfstash"]