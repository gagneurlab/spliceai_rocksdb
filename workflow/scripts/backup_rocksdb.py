import tarfile
import rocksdb


db = rocksdb.DB(snakemake.input['db'], 
                rocksdb.Options(max_open_files=5000), 
                read_only=True)
backup = rocksdb.BackupEngine(snakemake.output['backup'])
backup.create_backup(db, flush_before_backup=True)

tar = tarfile.open(snakemake.output['backup_gzip'], "w:gz")
tar.add(snakemake.output['backup'], arcname="backup")
tar.close()
