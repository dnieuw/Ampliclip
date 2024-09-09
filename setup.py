from setuptools import setup, find_packages
from ampliclip import __version__

setup(name='squirrel',
    version=__version__,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'ampliclip=ampliclip.__main__:main',
        ]
    },
    install_requires=[
        'biopython>=1.70',
        'pysam>=0.20'
    ],
    description='Tool to softclip reads in bam files based on amplicon primers',
    url='https://github.com/dnieuw/ampliclip',
    author='David F. Nieuwenhuijse',
    author_email='d.nieuwenhuijse@erasmusmc.nl',
    license='BSD 3-Clause',
    zip_safe=False)