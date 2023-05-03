from click.testing import CliRunner
from spliceai_rocksdb.spliceAI import SpliceAIDB
from spliceai_rocksdb.main import spliceai_rocksdb_download


def tests_download_spliceai_rocksdb(tmp_path, script_runner):
    runner = CliRunner()
    result = runner.invoke(spliceai_rocksdb_download,
                           f'--version _test --db_path {str(tmp_path)}')
    assert result.exit_code == 0

    db_download = SpliceAIDB(str(tmp_path))
    it = db_download.db.iterkeys()
    it.seek_to_first()
    download_variants = list(it)

    db = SpliceAIDB('tests/data/test_db')
    it = db.db.iterkeys()
    it.seek_to_first()
    variant = list(it)

    assert len(download_variants) == len(variant)
