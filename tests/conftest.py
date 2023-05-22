import pytest
from kipoiseq import Variant
from cyvcf2 import Writer
from kipoiseq.extractors import MultiSampleVCF


fasta_file = 'tests/data/hg19.nochr.chr17.fa'
multi_vcf_file = 'tests/data/multi_test.vcf.gz'
spliceai_snv_vcf_file_chr1 = 'tests/data/spliceai_snv_chr1.vcf'
spliceai_snv_vcf_file_chr9 = 'tests/data/spliceai_snv_chr9.vcf'
spliceai_rocksdb_chr1_prec = 'tests/data/spliceAI_chr1.db/'
spliceai_rocksdb_chr9_prec = 'tests/data/spliceAI_chr9.db/'


@pytest.fixture
def spliceai_snv_vcf(tmp_path):
    vcf_path = tmp_path / 'snv.vcf'

    header = '\n'.join([
        '##fileformat=VCFv4.2',
        '##contig=<ID=1,length=249250621>',
        '##contig=<ID=9,length=141213431>',
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">',

        '\t'.join(['#CHROM', 'POS', 'ID', 'REF',
                   'ALT', 'QUAL', 'FILTER', 'INFO'])
    ])

    writer = Writer.from_string(str(vcf_path), header)

    variants = [
        ('1:69091:A>C', 'SpliceAI=C|OR4F5|0.01|0.00|0.00|0.00|42|25|24|2'),
        ('1:69091:A>G', 'SpliceAI=G|OR4F5|0.00|0.00|0.07|0.00|43|42|-1|2'),
        ('9:37783955:A>C', 'SpliceAI=C|EXOSC3|0.00|0.00|0.00|0.00|0|-13|-44|-12'),
        ('9:37783955:A>C', 'SpliceAI=C|RP11-613M10.9|0.00|0.00|0.00|0.00|0|-13|-44|-12'),
        ('9:37783955:A>G', 'SpliceAI=G|EXOSC3|0.00|0.00|0.00|0.00|-11|-13|-44|24'),
        ('9:37783955:A>G', 'SpliceAI=G|RP11-613M10.9|0.00|0.00|0.00|0.00|38|-13|-44|24'),
        ('9:37783955:A>C', 'SpliceAI=C|test|0.00|0.00|0.00|0.00|0|-13|-44|-12')
    ]

    for v, info in variants:
        v = Variant.from_str(v)

        variant = writer.variant_from_string('\t'.join([
            v.chrom, str(v.pos), '.', v.ref, v.alt, '.', '.', info
        ]))

        writer.write_record(variant)
    writer.close()

    return vcf_path
