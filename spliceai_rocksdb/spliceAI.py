"""
This file provides python interface for spliceAI
"""
import rocksdb
from collections import namedtuple
from tqdm import tqdm
import numpy as np
import pandas as pd
from kipoiseq import Variant
from spliceai.utils import Annotator, get_delta_scores
from kipoiseq.extractors import MultiSampleVCF
import os
from pathlib import Path
import pathlib
from itertools import chain


def df_batch_writer(df_iter, output):
    df = next(df_iter)
    with open(output, 'w') as f:
        df.to_csv(f, index=False)

    for df in df_iter:
        with open(output, 'a') as f:
            df.to_csv(f, index=False, header=False)


def df_batch_writer_parquet(df_iter, output_dir, batch_size_parquet=100000):
    if not os.path.isdir(output_dir):
        output_dir.mkdir(exist_ok=True)

    dfs = list()
    num_rows = 0

    batch_num = 0
    for batch_num, df in enumerate(df_iter):
        dfs.append(df)
        num_rows += df.shape[0]
        if num_rows >= batch_size_parquet:
            df_all = pd.concat(dfs, axis=0)
            part_file = output_dir / f"{batch_num}.parquet"
            df_all.to_parquet(part_file, index=False, engine='pyarrow')
            dfs = list()
            num_rows = 0
    if num_rows > 0:
        batch_num += 1
        df_all = pd.concat(dfs, axis=0)
        part_file = output_dir / f"{batch_num}.parquet"
        df_all.to_parquet(part_file, index=False, engine='pyarrow')


def df_batch_writer_vcf(df_iter, output_path, header):
    import cyvcf2

    columns = [
        'alt_acceptor',
        'alt_acceptorIntron',
        'alt_donor',
        'alt_donorIntron',
        'alt_exon',
        'delta_logit_psi',
        'pathogenicity',
        'ref_acceptor',
        'ref_acceptorIntron',
        'ref_donor',
        'ref_donorIntron',
        'ref_exon'
    ]

    # vcf.add_info_to_header({
    #     'ID': 'SpliceAI',
    #     'Description': (
    #         "SpliceAIv1.3 variant annotation. "
    #         "These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), "
    #         "donor gain (DG), and donor loss (DL). "
    #         "Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"
    #     ),
    #     'Type': 'String',
    #     'Number': '.'
    # })

    vcf_writer = cyvcf2.Writer.from_string(str(output_path), header)
    for batch_df in df_iter:
        for idx, row in batch_df.iterrows():
            v = Variant.from_str(row.variant)
            vcf_variant = vcf_writer.variant_from_string('\t'.join([
                v.chrom, str(v.pos), '.', v.ref, v.alt, '.', '.', '.'
            ]))
            vcf_variant.INFO["SpliceAI"] = "|".join([
                v.alt, # ALLELE
                row["gene_name"], # SYMBOL
                str(row["acceptor_gain"]), # DS_AG
                str(row["acceptor_loss"]), # DS_AL
                str(row["donor_gain"]), # DS_DG
                str(row["donor_loss"]), # DS_DL
                str(row["acceptor_gain_position"]), # DP_AG
                str(row["acceptor_loss_position"]), # DP_AL
                str(row["donor_gain_position"]), # DP_DG
                str(row["donor_loss_position"]), # DP_DL
            ])
            vcf_writer.write_record(vcf_variant)
    vcf_writer.close()


class VariantDB:

    def __init__(self, path):
        self.db = rocksdb.DB(
            path,
            rocksdb.Options(
                create_if_missing=True,
                max_open_files=300,
            ),
            read_only=True
        )


    @staticmethod
    def _variant_to_byte(variant):
        return bytes(str(variant), 'utf-8')

    def _type(self, value):
        raise NotImplementedError()

    def _get(self, variant):
        if variant.startswith('chr'):
            variant = variant[3:]
        return self.db.get(self._variant_to_byte(variant))

    def __getitem__(self, variant):
        value = self._get(variant)
        if value:
            return self._type(value)
        else:
            raise KeyError('This variant "%s" is not in the db'
                           % str(variant))

    def __contains__(self, variant):
        return self._get(variant) is not None

    def get(self, variant, default=None):
        try:
            return self[variant]
        except KeyError:
            return default

    def items(self):
        it = self.db.iteritems()
        it.seek_to_first()
        for variant, value in it:
            yield variant.decode('utf-8'), self._type(value)


class SpliceAIDB(VariantDB):

    def _type(self, value):
        return list(self._parse(value))

    def _parse(self, value):
        for i in value.decode('utf-8').split(';'):
            results = i.split('|')
            scores = np.array(list(map(float, results[1:])))
            yield SpliceAI.Score(
                results[0], scores[:4].max(), *scores
            )


class SpliceAI:
    Score = namedtuple('Score', ['gene_name', 'delta_score',
                                 'acceptor_gain', 'acceptor_loss',
                                 'donor_gain', 'donor_loss',
                                 'acceptor_gain_position',
                                 'acceptor_loss_position',
                                 'donor_gain_position',
                                 'donor_loss_position'])
    Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])

    def __init__(self, fasta=None, annotation=None, db_path=None,
                 dist=50, mask=1):
        """
        Args:
          fasta: fasta file path
          annotation: 'grch37' or 'grch38'
          dist: area of interest based on distance to variant
          mask: mask for 'N'
        """
        assert ((fasta is not None) and (annotation is not None)) \
               or (db_path is not None)
        self.db_only = fasta is None
        self.annotation = str(annotation).lower()
        if self.annotation not in {"grch37", "grch38"}:
            raise ValueError(f"Unknown annotation version: '{annotation}'!")
        if not self.db_only:
            self.annotator = Annotator(fasta, annotation)
        self.dist = dist
        self.mask = mask
        self.db = {} if db_path else None
        if db_path:
            for chr in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']:
                try:
                    self.db[chr] = SpliceAIDB(db_path[chr])
                except KeyError:
                    print(f"There was no database for chr{chr} provided")
                    continue

    @staticmethod
    def _to_record(variant):
        if type(variant) == str:
            variant = Variant.from_str(variant)
        return SpliceAI.Record(variant.chrom, variant.pos,
                               variant.ref, [variant.alt])

    @staticmethod
    def parse(output):
        results = output.split('|')
        results = [0 if i == '.' else i for i in results]
        scores = np.array(list(map(float, results[2:])))
        return SpliceAI.Score(
            results[1], scores[:4].max(), *scores
        )

    def predict(self, variant):
        record = self._to_record(variant)
        if self.db:
            try:
                if record[0].startswith('chr'):
                    return self.db[record[0][3:]][str(variant)]
                return self.db[record[0]][str(variant)]
            except KeyError:
                if self.db_only:
                    return []
        return [
            self.parse(i)
            for i in get_delta_scores(record, self.annotator,
                                      self.dist, self.mask)
        ]

    def predict_df(self, variants, vcf=None):
        rows = self._predict_df(variants, vcf)
        columns = [
            'variant', 'gene_name', 'delta_score',
            'acceptor_gain', 'acceptor_loss',
            'donor_gain', 'donor_loss',
            'acceptor_gain_position',
            'acceptor_loss_position',
            'donor_gain_position',
            'donor_loss_position'
        ]
        type_dict = {
            'variant': 'string',
            'gene_name': 'string',
            'delta_score': 'float64',
            'acceptor_gain': 'float64',
            'acceptor_loss': 'float64',
            'donor_gain': 'float64',
            'donor_loss': 'float64',
            'acceptor_gain_position': 'int64',
            'acceptor_loss_position': 'int64',
            'donor_gain_position': 'int64',
            'donor_loss_position': 'int64',
        }
        return pd.DataFrame(rows, columns=columns).astype(type_dict).set_index('variant')

    def _predict_df(self, variants, vcf=None):
        for v in variants:
            for score in self.predict(v):
                row = score._asdict()
                row['variant'] = str(v)
                yield row

    def _predict_on_vcf(self, vcf_file, batch_size=100000):
        vcf = MultiSampleVCF(vcf_file)
        for variants in vcf.batch_iter(batch_size):
            yield self.predict_df(variants, vcf).reset_index('variant')

    def predict_save(self, vcf_file, output_path,
                     batch_size=100000, progress=True):
        batches = self._predict_on_vcf(vcf_file, batch_size=batch_size)

        if progress:
            batches = iter(tqdm(batches))

        if not isinstance(output_path, pathlib.PosixPath):
            output_path = Path(output_path)

        file_ext = output_path.suffix.lower()

        if file_ext == '.csv':
            df_batch_writer(batches, output_path)
        elif file_ext == '.parquet':
            df_batch_writer_parquet(batches, output_path, batch_size)
        elif file_ext == '.vcf':
            from spliceai_rocksdb import header

            if self.annotation == "grch37":
                vcf_header = header.header_grch37
            elif self.annotation == "grch38":
                vcf_header = header.header_grch38
            else:
                # this part of the code should not be reachable
                raise ValueError(f"Unknown annotation version: '{self.annotation}'!")

            df_batch_writer_vcf(batches, output_path, header=vcf_header)
        else:
            raise ValueError(
                'File extension %s is not supported. '
                'Supported file types are: csv, vcf, parquet' % file_ext)
