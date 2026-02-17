from setuptools import setup

setup(name='ViralFlow',
      version='1.3.6',
      description='''
      Nextflow workflow for reference-based viral genome assembly, quality control, 
      variant calling, and lineage assignment. Supports multiple viruses with 
      containerized tools for reproducible genomic surveillance.
      ''',
      url='https://github.com/dezordi/ViralFlow/',
      author='Antonio Marinho & Filipe Z. Dezordi',
      author_email='amarinhosn@pm.me & zimmer.filipe@gmail.com',
      py_modules = [],
      classifiers=[
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
          ],
      scripts=['wrapper/viralflow'],
      zip_safe=False
)
