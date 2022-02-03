import pytest
import pandas as pd
import rocksdb
from spliceai_rocksdb.spliceAI import SpliceAI
from conftest import fasta_file, multi_vcf_file


@pytest.fixture
def spliceai():
    return SpliceAI(fasta_file, 'grch37')


def test_SpliceAI_parse():
    score = SpliceAI.parse('CA|TTN|0.07|1.00|0.00|0.00|-7|-1|35|-29')

    assert score == SpliceAI.Score(
        'TTN', 1., 0.07, 1., 0., 0., -7, -1, 35, -29)


def test_SpliceAI_predict(spliceai):
    scores = spliceai.predict('17:34149615:A>T')
    # 'TAF15|0.25|0|0|0|11|29|11|33'

    assert scores[0] == spliceai.Score(
        gene_name='TAF15', delta_score=0.25,
        acceptor_gain=0.25, acceptor_loss=0.0,
        donor_gain=0.0, donor_loss=0.0,
        acceptor_gain_position=11, acceptor_loss_position=29,
        donor_gain_position=11, donor_loss_position=33)

    scores = spliceai.predict('17:883700:A>T')
    assert scores == []


def test_SpliceAI_predict_df(spliceai):
    df = spliceai.predict_df([
        '17:34149615:A>T',
        '17:76202989:CAGATTATCTTGAA>C',
        '17:883700:A>T'
    ])
    pd.testing.assert_frame_equal(df, pd.DataFrame({
        'variant': pd.Series(['17:34149615:A>T', '17:76202989:CAGATTATCTTGAA>C'], dtype='string'),
        'gene_name': pd.Series(['TAF15', 'AFMID'], dtype='string'),
        'delta_score': [0.25, 0.98],
        'acceptor_gain': [0.25, 0.74],
        'acceptor_loss': [0.0, 0.98],
        'donor_gain': [0.0, 0.0],
        'donor_loss': [0.0, 0.0],
        'acceptor_gain_position': [11, 27],
        'acceptor_loss_position': [29, 3],
        'donor_gain_position': [11, 30],
        'donor_loss_position': [33, 3]
    }).set_index('variant'))


@pytest.fixture
def spliceai_rocksdb(tmp_path):
    db_path = str(tmp_path / 'db')
    db = rocksdb.DB(db_path,
                    rocksdb.Options(create_if_missing=True))

    data = {
        '17:34149615:A>T': 'X|0.07|1.00|0.00|0.00|-7|-1|35|-29',
    }

    for k, v in data.items():
        db.put(bytes(k, 'utf-8'), bytes(v, 'utf-8'))

    return db_path


@pytest.fixture
def spliceai_db(spliceai_rocksdb):
    return SpliceAI(fasta_file, 'grch37', db_path=spliceai_rocksdb)


def test_SpliceAI_predict_with_db(spliceai_db):
    scores = spliceai_db.predict('17:34149615:A>T')
    assert scores[0] == spliceai_db.Score(
        gene_name='X', delta_score=1.0,
        acceptor_gain=0.07, acceptor_loss=1.0,
        donor_gain=0.0, donor_loss=0.0,
        acceptor_gain_position=-7, acceptor_loss_position=-1,
        donor_gain_position=35, donor_loss_position=-29)


def test_SpliceAI_predict_on_vcf(spliceai_db, tmp_path):
    output_csv = tmp_path / 'output.csv'
    spliceai_db.predict_save(multi_vcf_file, output_csv)
    df = pd.read_csv(output_csv)

    assert df.columns.tolist() == [
        'variant', 'gene_name', 'delta_score', 'acceptor_gain',
        'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
        'acceptor_loss_position', 'donor_gain_position', 'donor_loss_position'
    ]

    output_parquet = tmp_path / 'output.parquet'
    spliceai_db.predict_save(multi_vcf_file, output_parquet)
    df = pd.read_parquet(output_parquet, engine='pyarrow')

    assert df.columns.tolist() == [
        'variant', 'gene_name', 'delta_score', 'acceptor_gain',
        'acceptor_loss', 'donor_gain', 'donor_loss', 'acceptor_gain_position',
        'acceptor_loss_position', 'donor_gain_position', 'donor_loss_position'
    ]

    output_vcf = tmp_path / 'output.vcf'
    spliceai_db.predict_save(multi_vcf_file, output_vcf)
    import cyvcf2
    vcf = cyvcf2.VCF(output_vcf)
    variant = next(vcf)
    assert 'CA|BRCA1|0.0|0.0|0.0|0.0|0|0|0|0' == variant.INFO.get("SpliceAI")
    vcf.close()




def test_SpliceAI_predict_db_only(spliceai_rocksdb):
    spliceai = SpliceAI(db_path=spliceai_rocksdb)
    df = spliceai.predict_df(['17:34149615:A>T', '17:34149617:A>T'])
    assert df.shape == (1, 10)

    spliceai = SpliceAI(db_path=spliceai_rocksdb)
    df = spliceai.predict_df(['chr17:34149615:A>T', 'chr17:34149617:A>T'])
    assert df.shape == (1, 10)
