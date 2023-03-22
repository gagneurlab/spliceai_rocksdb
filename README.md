# SpliceAI Rocksdb

Fast look up interface and prediction for missing variants in look table for spliceAI with rocksdb.

## Installation

Clone this git repository:

```
git clone https://github.com/gagneurlab/spliceai_rocksdb.git
```

Go to the cloned repository:

```
cd spliceai_rocksdb
```

Create conda environment:

```
# Recommended if you have mamba installed
mamba env create -f environment.yml
# otherwise
conda env create -f environment.yml
```

Activate conda environment:

```
conda activate spliceai-rocksdb
```

Install modules from spliceai_rocksdb

```
pip install -e .
```

## Download database

Download rocksdb for gnomad
```
spliceai_rocksdb_download --version {version} --db_path {output_path}
```

Supported version (grch37, grch38)


## Usage

To run spliceai_rocksdb in the command line use:
```
spliceai_rocksdb -I input.vcf -O output.[csv, parquet, vcf] -R genome.fa -A grch37 -db {db_path}
```

Alternatively use the snakemake pipeline in /workflow 
```
cd workflow
# modify config.yaml and setup datasets
python -m snakemake -j 1
```