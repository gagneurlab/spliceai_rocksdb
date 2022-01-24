import rocksdb
from tqdm import tqdm
from kipoiseq.extractors import MultiSampleVCF


def create_spliceai_rocksdb(db_path, vcfs, max_open_files=300,
                            create_if_missing=True, batch_size=10000):
    options = rocksdb.Options(
        create_if_missing=create_if_missing,
        max_open_files=max_open_files)
    db = rocksdb.DB(db_path, options)

    for vcf_path in vcfs:
        vcf = MultiSampleVCF(vcf_path)

        for batch in tqdm(vcf.batch_iter(batch_size)):
            batch_variants = dict()

            for variant in batch:
                if str(variant) not in batch_variants:
                    _row = db.get(bytes(str(variant), 'utf-8'))
                    if _row:
                        batch_variants[str(variant)] = [_row.decode("utf-8")]
                    else:
                        batch_variants[str(variant)] = list()

                row = variant.source.INFO.get('SpliceAI')
                row = '|'.join(row.split('|')[1:])
                batch_variants[str(variant)].append(row)

            batch_writer = rocksdb.WriteBatch()
            for variant, rows in batch_variants.items():
                variant = bytes(variant, 'utf-8')
                batch_writer.put(variant, bytes(';'.join(rows), 'utf-8'))
            db.write(batch_writer)
