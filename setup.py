import os
from setuptools import setup, find_packages


# read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()


setup(
    name='metalink',
    version='0.3.7',
    url='https://github.com/LiuyangLee/Metalink',
    license='GPL-3',
    author='Liuyang Li',
    author_email='cyanobacteria@yeah.net',
    description=('Metalink: Joint taxonomic and functional profiling of metagenomes and metatranscriptomes across resolution scales. Metalink can be run in different modes to obtain joint outputs of functional and taxonomic profiles at multiple levels'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['metalink'],
    include_package_data=True,  # include all files in MANIFEST.in
    data_files=[],
    zip_safe=False,
    install_requires=[
        'click>=8.0,<9.0',  # uses ctx.get_parameter_source / click.core.ParameterSource
        'ruamel.yaml',       # used for config loading/dumping
        'snakemake>=7.32.4',
    ],
    # Note: most heavy bioinformatics dependencies are installed via conda.
    # The click pin also protects the use of click.formatting helpers and
    # __click_params__ internals in metalink.py.
    entry_points={
          'console_scripts': [
              'metalink = metalink.metalink:cli'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
