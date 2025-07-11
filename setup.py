from setuptools import setup

setup(name='ViralFlow',
      version='1.3.1',
      description='''
      Workflows for viral genome Assembly at FioCruz/IAM
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
