from spliceai_rocksdb.spliceAI import SpliceAI


if snakemake.params['lookup_only']:
    model = SpliceAI(db_path=snakemake.input['db'])
else:
    model = SpliceAI(snakemake.input['fasta'],
                     annotation=snakemake.params['genome'],
                     db_path=snakemake.input['db'])


model.predict_save(snakemake.input['vcf'],
                   snakemake.output['csv'])
