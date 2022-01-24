from spliceai_rocksdb.create import create_spliceai_rocksdb


vcfs = [
    snakemake.input['snv'],
    snakemake.input['indel']
]
create_spliceai_rocksdb(snakemake.output['db'], vcfs,
                        batch_size=snakemake.params['batch_size'])
