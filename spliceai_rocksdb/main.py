import rocksdb
import wget
import click
import tarfile
from spliceai_rocksdb import SpliceAI

db_url = {
    'grch37': 'https://nextcloud.in.tum.de/index.php/s/dEDnN8gBEsAbtCq/download',
    'grch38': '',
    '_test': ''
}


@click.command()
@click.option('--version', help='SpliceAI rocksdb version (currently grch37, grch38 supported)')
@click.option('--db_path', help='Path to download database')
def spliceai_rocksdb_download(version, db_path):

    if version not in db_url:
        raise(f'Version {version} is not supported.')

    print('Downloading database...')
    download_path = db_path + '_backup.tar.gz'
    wget.download(db_url[version], out=download_path)

    print('\nUnzipping database...')
    file = tarfile.open(download_path)
    backup_dir = db_path + '_backup/'
    file.extractall(backup_dir)
    file.close()

    print('Storing database...')
    backup = rocksdb.BackupEngine(backup_dir + 'backup/')
    backup.restore_latest_backup(db_path, db_path)


@click.command()
@click.option('-A', help='Genome annotation: `grch37`, `grch38`')
@click.option('-R', help='Path for fasta file `grch37`, `grch38`')
@click.option('--I', help='Path for the input VCF')
@click.option('--db', help='Path for the spliceAI rocksdb')
@click.option('--I', help='Output file path')
@click.option('--look_up_only', help='Performs only variant look up.'
              'Do not performs if variants are missing in database',
              is_flag=True)
@click.option('--batch_size', help='Batch size while performing predictions',
              default=100000, show_default=True)
def spliceai_rocksdb(annotation, fasta_file, vcf, db, output,
                     look_up_only, batch_size):
    if look_up_only:
        spliceai = SpliceAI(None, annotation, db)
    spliceai.predict_save(vcf, output, batch_size=batch_size)
