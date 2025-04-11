from setuptools import setup, find_packages

setup(
    name="vcfstash",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "vcfstash=src.cli:main",
        ],
    },
    install_requires=[
        "pysam",
        "pytest"
    ],
    extras_require={
        "parquet": ["pandas", "pyarrow"],
    },
)