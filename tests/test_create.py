import rocksdb
from spliceai_rocksdb.create import create_spliceai_rocksdb


def test_create_spliceai_rocksdb(tmp_path, spliceai_snv_vcf):
    db_path = str(tmp_path / 'db')

    create_spliceai_rocksdb(db_path, [spliceai_snv_vcf], batch_size=2)

    db = rocksdb.DB(
        db_path,
        rocksdb.Options(
            create_if_missing=True,
            max_open_files=300
        ),
        read_only=True
    )

    variant = bytes(str('9:37783955:A>C'), 'utf-8')
    value = db.get(variant).decode('utf-8').split(';')

    assert len(value) == 3
