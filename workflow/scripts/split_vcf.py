import os
from pathlib import Path
from kipoiseq.extractors import MultiSampleVCF

vcf = MultiSampleVCF(snakemake.input['vcf'])

batch_size = snakemake.params['batch_size']

for batch in vcf.query_all().batch_iter(batch_size):
    _, interval = batch.variant_intervals[0]
    interval = f'{interval.chrom}:{interval.start}-{interval.end}'
    prefix = snakemake.params['prefix_vcf']
    path = Path(snakemake.params['vcfs'])
    path.mkdir(exist_ok=True)
    path = path / f'{prefix}_{interval}.vcf.gz'
    batch.to_vcf(path, remove_samples=True, clean_info=True)
