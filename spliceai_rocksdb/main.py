import rocksdb
import pooch
import click
import tarfile
from spliceai_rocksdb import SpliceAI

# db_url = {
#     'grch37': 'https://nextcloud.in.tum.de/index.php/s/dEDnN8gBEsAbtCq/download',
#     'grch38': 'https://nextcloud.in.tum.de/index.php/s/SzHjFkd9Y6jp3Dk/download',
#     '_test': ''
# }

db_url = {
    'grch37': {
        '1': ['doi:10.5281/zenodo.7925611/spliceAI_rocksdb_hg19_chr1.tar.gz', 'md5:a70189ed3da20e76ce61540c77263035'],
		'2': ['doi:10.5281/zenodo.7925735/spliceAI_rocksdb_hg19_chr2.tar.gz', 'md5:7db7a23ae197a5c86cc65cce4d9799b1'],
		'3': ['doi:10.5281/zenodo.7925768/spliceAI_rocksdb_hg19_chr3.tar.gz', 'md5:f661c1fe1301d1fe5ef6da7a53864a24'],
		'4': ['doi:10.5281/zenodo.7925900/spliceAI_rocksdb_hg19_chr4.tar.gz', 'md5:8b4eb2b036a3ced76001249eef21335e'],
		'5': ['doi:10.5281/zenodo.7925902/spliceAI_rocksdb_hg19_chr5.tar.gz', 'md5:790ed3735ebe90ffa548a12a3412b99e'],
		'6': ['doi:10.5281/zenodo.7925908/spliceAI_rocksdb_hg19_chr6.tar.gz', 'md5:25746fd6a62eb13b39c57c8b9de666c7'],
		'7': ['doi:10.5281/zenodo.7925923/spliceAI_rocksdb_hg19_chr7.tar.gz', 'md5:b1e57191d887dc6ffc74db28685aad6c'],
		'8': ['doi:10.5281/zenodo.7925931/spliceAI_rocksdb_hg19_chr8.tar.gz', 'md5:33cf33ed6a0b930009400ccace0af4f9'],
		'9': ['doi:10.5281/zenodo.7925949/spliceAI_rocksdb_hg19_chr9.tar.gz', 'md5:42315cf123009d363d2cb6c97e07508c'],
		'10': ['doi:10.5281/zenodo.7925959/spliceAI_rocksdb_hg19_chr10.tar.gz', 'md5:9d2fb68b40e90452686c8cdbfa22d0b4'],
		'11': ['doi:10.5281/zenodo.7925967/spliceAI_rocksdb_hg19_chr11.tar.gz', 'md5:663057ce569a5270b217cae56052b748'],
		'12': ['doi:10.5281/zenodo.7925977/spliceAI_rocksdb_hg19_chr12.tar.gz', 'md5:40c656174ceeed5f30d6e37cce170274'],
		'13': ['doi:10.5281/zenodo.7925984/spliceAI_rocksdb_hg19_chr13.tar.gz', 'md5:1dd93a6458ad55ebbc7a30c3436f0e27'],
		'14': ['doi:10.5281/zenodo.7925993/spliceAI_rocksdb_hg19_chr14.tar.gz', 'md5:f20e62af6ea7787c651435f49d13eea9'],
		'15': ['doi:10.5281/zenodo.7926008/spliceAI_rocksdb_hg19_chr15.tar.gz', 'md5:63b68f1b8050a733c0dacfe91b18432f'],
		'16': ['doi:10.5281/zenodo.7926021/spliceAI_rocksdb_hg19_chr16.tar.gz', 'md5:af001d5522822b717e78a80a098c6838'],
		'17': ['doi:10.5281/zenodo.7926028/spliceAI_rocksdb_hg19_chr17.tar.gz', 'md5:7c3d9698b09bad5fa2c87b6a622ee54b'],
		'18': ['doi:10.5281/zenodo.7926032/spliceAI_rocksdb_hg19_chr18.tar.gz', 'md5:5b2d01d90a97eb2771f4dae75e598e45'],
		'19': ['doi:10.5281/zenodo.7926040/spliceAI_rocksdb_hg19_chr19.tar.gz', 'md5:8f50f01b5b34e29a36ca07809e9462f6'],
		'20': ['doi:10.5281/zenodo.7926052/spliceAI_rocksdb_hg19_chr20.tar.gz', 'md5:a39d04b257d8134a6e8a39cbd2bcd844'],
		'21': ['doi:10.5281/zenodo.7926058/spliceAI_rocksdb_hg19_chr21.tar.gz', 'md5:cbd5b24faa37847f08f4274795ae999a'],
		'22': ['doi:10.5281/zenodo.7926064/spliceAI_rocksdb_hg19_chr22.tar.gz', 'md5:fedbf9bb7509b5bc4d97becd390f08a6'],
		'X': ['doi:10.5281/zenodo.7926068/spliceAI_rocksdb_hg19_chrX.tar.gz', 'md5:9eb5ff37219bc1caf03a38b4b28dd247'],
		'Y': ['doi:10.5281/zenodo.7926074/spliceAI_rocksdb_hg19_chrY.tar.gz', 'md5:55b18fb7eed154591a3d4d1ee2f24248']
        },
    'grch38': {
        '1': ['doi:10.5281/zenodo.7926110/spliceAI_rocksdb_hg38_chr1.tar.gz', 'md5:75bbd5b1d203c0e322daf071df4b5a73'],
 		'2': ['doi:10.5281/zenodo.7926108/spliceAI_rocksdb_hg38_chr2.tar.gz', 'md5:f41868b7cf95dd3d1b829906cf4e32f8'],
 		'3': ['doi:10.5281/zenodo.7926124/spliceAI_rocksdb_hg38_chr3.tar.gz', 'md5:bc591cc688f785806c8fe9559c99fd06'],
 		'4': ['doi:10.5281/zenodo.7926133/spliceAI_rocksdb_hg38_chr4.tar.gz', 'md5:7e8e7e0b92b38fc654fdcdb11769e449'],
 		'5': ['doi:10.5281/zenodo.7926126/spliceAI_rocksdb_hg38_chr5.tar.gz', 'md5:ca142b0a942f55749831894cc77ad278'],
 		'6': ['doi:10.5281/zenodo.7926137/spliceAI_rocksdb_hg38_chr6.tar.gz', 'md5:16f30f547071cdceb37800234b06fdf8'],
 		'7': ['doi:10.5281/zenodo.7926141/spliceAI_rocksdb_hg38_chr7.tar.gz', 'md5:3138c7e65ba340cd4b68341e81ac8998'],
 		'8': ['doi:10.5281/zenodo.7926145/spliceAI_rocksdb_hg38_chr8.tar.gz', 'md5:d4e72bf6a3813a3d22e953ba516360cb'],
 		'9': ['doi:10.5281/zenodo.7926149/spliceAI_rocksdb_hg38_chr9.tar.gz', 'md5:a76826992511b99bf4bf6f4b3ac9c218'],
 		'10': ['doi:10.5281/zenodo.7926157/spliceAI_rocksdb_hg38_chr10.tar.gz', 'md5:26fd6f35e11f82a3c5d2be56e6c1f9cd'],
 		'11': ['doi:10.5281/zenodo.7928451/spliceAI_rocksdb_hg38_chr11.tar.gz', 'md5:8a9910858c1e7daad67a8e52a810800a'],
 		'12': ['doi:10.5281/zenodo.7928469/spliceAI_rocksdb_hg38_chr12.tar.gz', 'md5:3cac99dccfa3cac8fecbf0474425a095'],
 		'13': ['doi:10.5281/zenodo.7928475/spliceAI_rocksdb_hg38_chr13.tar.gz', 'md5:228f0d714e371bf97cf76913270bbd7c'],
 		'14': ['doi:10.5281/zenodo.7928489/spliceAI_rocksdb_hg38_chr14.tar.gz', 'md5:9643188534b6072dec17e8d243ba554e'],
 		'15': ['doi:10.5281/zenodo.7928480/spliceAI_rocksdb_hg38_chr15.tar.gz', 'md5:88f0c938ec93b0b0593b3560caf45168'],
 		'16': ['doi:10.5281/zenodo.7928478/spliceAI_rocksdb_hg38_chr16.tar.gz', 'md5:3a734b9d1bfadf514db73f1316ee4f6b'],
 		'17': ['doi:10.5281/zenodo.7928518/spliceAI_rocksdb_hg38_chr17.tar.gz', 'md5:7db176648c4999a6c2152404c5bf2681'],
 		'18': ['doi:10.5281/zenodo.7928526/spliceAI_rocksdb_hg38_chr18.tar.gz', 'md5:a7b87ebb77adc5d1fef782e5e18ac328'],
 		'19': ['doi:10.5281/zenodo.7928532/spliceAI_rocksdb_hg38_chr19.tar.gz', 'md5:065fe778278fd7d84790d65b3160d65b'],
 		'20': ['doi:10.5281/zenodo.7928537/spliceAI_rocksdb_hg38_chr20.tar.gz', 'md5:e34e2caa19a9ab79adad0cef86762af2'],
 		'21': ['doi:10.5281/zenodo.7928539/spliceAI_rocksdb_hg38_chr21.tar.gz', 'md5:b625796d6c501c25c0c5c1bcdf37c539'],
 		'22': ['doi:10.5281/zenodo.7928543/spliceAI_rocksdb_hg38_chr22.tar.gz', 'md5:b5b9640b401942e47a9e2a4e78c60914'],
 		'X': ['doi:10.5281/zenodo.7928552/spliceAI_rocksdb_hg38_chrX.tar.gz', 'md5:a915558ca341f44c1b8ad44198423e9b'],
 		'Y': ['doi:10.5281/zenodo.7928556/spliceAI_rocksdb_hg38_chrY.tar.gz', 'md5:c3e8b0404dad1abdeeb259e2288595ec']
        },
    '_test': {
        '1': ['doi:10.5281/zenodo.7990703/SpliceAI_rocksdb_hg19_test_chr1.tar.gz', 'md5:803e44d671d89eb17590e85524719669'],
        '2': ['doi:10.5281/zenodo.7990703/SpliceAI_rocksdb_hg19_test_chr2.tar.gz', 'md5:a2fedb977c7d7ee1747e4c408027044c']
	}
}

@click.command()
@click.option('--version', help='SpliceAI rocksdb version (currently grch37, grch38 supported)')
@click.option('--db_path', help='Path to download database')
@click.option('--chromosome', help='Number of chromosome for which database will be downloaded')
def spliceai_rocksdb_download(version, db_path, chromosome):
    
    if chromosome.startswith('chr'):
        chromosome = chromosome[3:]

    if version not in db_url:
        raise(f'Version {version} is not supported.')

    print('Downloading database...')
    download_path = db_path + '_backup.tar.gz'

    filename = pooch.retrieve(
        url=db_url[version][chromosome][0],
        known_hash=db_url[version][chromosome][1],
        fname=download_path.split('/')[-1],
        path='/'.join(download_path.split('/')[:-1]) + '/'
        )

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
