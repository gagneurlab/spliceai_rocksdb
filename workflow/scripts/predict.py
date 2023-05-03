from spliceai_rocksdb.spliceAI import SpliceAI


if snakemake.params['lookup_only']:
    model = SpliceAI(db_path='data/results/dbs_by_chrom/spliceAI_hg19_chr1.db/')
else:
    model = SpliceAI(snakemake.input['fasta'],
                     annotation=snakemake.params['genome'],
                     db_path='data/results/dbs_by_chrom/spliceAI_hg19_chr1.db/')


model.predict_save(snakemake.input['vcf'],
                   snakemake.output['csv'])
