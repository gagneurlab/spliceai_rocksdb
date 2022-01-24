from setuptools import setup, find_packages


requirements = [
    'setuptools',
    'click',
    'python-rocksdb',
    'wget',
    'spliceai',
    'kipoiseq',
    'cyvcf2',
    'tqdm'
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'pytest-console-scripts']

setup(
    author="M. Hasan Çelik",
    author_email='muhammedhasancelik@gmail.com',
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    description="Fast look up interface and prediction for SpliceAI with rocksdb.",
    install_requires=requirements,
    license="MIT license",
    entry_points={
        'console_scripts': [
            'spliceai_rocksdb_download=spliceai_rocksdb.main:spliceai_rocksdb_download',
            'spliceai_rocksdb=spliceai_rocksdb.main:spliceai_rocksdb'
        ]
    },
    keywords=['genomics', 'gnomad', 'variant',
              'spliceai', 'variant effect predict'],
    name='spliceai_rocksdb',
    packages=find_packages(include=['spliceai_rocksdb']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gagneurlab/spliceai_rocksdb',
    version='0.0.1'
)
