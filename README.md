# SpliceAI Rocksdb

Fast look up interface and prediction for missing variants in look table for spliceAI with rocksdb.

## Installation

```
conda install -c conda-forge rocksdb python-rocksdb
pip install spliceai_rocksdb
```

## Download database

Download rocksdb for gnomad
```
spliceai_rocksdb_download --version {version} --db_path {output_path}
```

Supported version (grch37, grch38)

## Run

```
spliceai_rocksdb -I input.[csv, parquet, vcf] -O output.vcf -R genome.fa -A grch37 -db {db_path}
```

## Create Database

```
pip install snakemake cython
# modify workflow/config.yaml and setup datasets
python -m snakemake -j 1
```