from distutils.core import setup
import glob

# # load version directly from the module
# version = __import__('genial').version
version = '0.1'

scripts = glob.glob('bin/*')
requirements = open('requirements.txt').readlines()

setup(
    name='genial',
    author='Bruno F Souza',
    author_email='fsouza.bruno@gmail.com',
    version=version,
    packages=['genial', 'genial/gff'],
    scripts=scripts,
    url='https://github.com/varnion/genial',
    license='BSD',
    install_requires=requirements,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]

      )
