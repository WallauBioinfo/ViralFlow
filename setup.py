from setuptools import setup, find_packages

setup(
    name='ViralFlow',
    version='1.4.0',
    description='''
    Workflows for viral genome Assembly at FioCruz/IAM
    ''',
    url='https://github.com/WallauBioinfo/ViralFlow',
    author='Antonio Marinho & Filipe Z. Dezordi',
    author_email='amarinhosn@pm.me & zimmer.filipe@gmail.com',
    packages=find_packages(),
    install_requires=[
        'click>=8.0',
    ],
    entry_points={
        'console_scripts': [
            'viralflow=wrapper.cli:cli',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    zip_safe=False
)
