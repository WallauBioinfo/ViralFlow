from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install
import subprocess
import sys
import os

class PostDevelopCommand(develop):
    """Command executed for 'pip install -e'"""
    def run(self):
        develop.run(self)
        print("\n" + "="*60)
        print("Configurating Singularity for ViralFlow...")
        print("="*60)
        # Get the directory of the config file
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        config_script = os.path.join(setup_dir, 'wrapper', 'configure_singularity.py')
        subprocess.check_call([sys.executable, config_script])

class PostInstallCommand(install):
    """Command executed for 'pip install'"""
    def run(self):
        install.run(self)
        print("\n" + "="*60)
        print("Configurating Singularity for ViralFlow...")
        print("="*60)
        # Get the directory of the setup.py   
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        config_script = os.path.join(setup_dir, 'wrapper', 'configure_singularity.py')
        subprocess.check_call([sys.executable, config_script])

setup(name='ViralFlow',
      version='1.3.5',
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
      cmdclass={
          'develop': PostDevelopCommand,
          'install': PostInstallCommand,
      },
      zip_safe=False
)