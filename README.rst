genial: GENome Interactive Annotation Library
=============================================

.. image:: https://badge.fury.io/py/genial.png
    :target: https://badge.fury.io/py/genial

This library provides `InteractiveAnnotation`, a high level representation of genome annotations,
allowing users to easily extract information, manipulate and reformat commonly used annotation
files such as BED and GFF.

The package is currently in alpha stage and only runs on `python3`.

Supported formats
-----------------

* Input: BED, GFF3, GTF

* Output: BED


Scripts
-------
For convenience, we provide two CLI utilities: `annotParser.py` and `annotMergeSmallGap.py`.

.. code-block:: bash

    $ annotParser.py -h
    usage: annotParser.py [-h] [-i INPUT] [-o OUTPUT] [-f {gff3,gtf,bed}]
                          [-t {extb,bed}] [-n MIN_EXON_COUNT]
                          [-igs IGNORE_GAPS_SMALLER_THAN]
                          [-igb IGNORE_GAPS_BIGGER_THAN] [-v]

    Parse, filter and convert annotation files

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            input file. to read from pipe, use the argument
                            'stdin'
      -o OUTPUT, --output OUTPUT
                            output file
      -f {gff3,gtf,bed}, --input_format {gff3,gtf,bed}
                            input file format
      -t {extb,bed}, --output_format {extb,bed}
                            output file format
      -n MIN_EXON_COUNT, --min_exon_count MIN_EXON_COUNT
                            min number of exons
      -igs IGNORE_GAPS_SMALLER_THAN, --ignore_gaps_smaller_than IGNORE_GAPS_SMALLER_THAN
      -igb IGNORE_GAPS_BIGGER_THAN, --ignore_gaps_bigger_than IGNORE_GAPS_BIGGER_THAN
      -v, --invert_match    select non matching annotations (similar to grep -v)



.. code-block:: bash

    $ annotMergeSmallGaps.py -h
    usage: annotMergeSmallGaps.py [-h] [-i INPUT] [-o OUTPUT] [-f {gff3,bed,gtf}]
                                  [-t {extb,bed}] [-s SMALL_GAP_SIZE]

    Merge exons separated by small gaps. Can also be used to convert different
    kinds of annotations.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            input file. to read from pipe, use the argument
                            'stdin'
      -o OUTPUT, --output OUTPUT
                            output file
      -f {gff3,bed,gtf}, --input_format {gff3,bed,gtf}
                            input file format
      -t {extb,bed}, --output_format {extb,bed}
                            output file format
      -s SMALL_GAP_SIZE, --small_gap_size SMALL_GAP_SIZE
                            gap size.



Both the scripts above can also be used to convert from different kinds of annotation files.
A more advanced usage can be achieved importing the library.

Instalation instructions
------------------------

for a systemwide instalation on ubuntu:

.. code-block:: bash

    sudo apt-get install python3-pip
    sudo pip3 install genial


Acknowledgements
----------------

I'd like to thank my friends, `Lucas Silva`_ and David Pires, for all the help and encouragement to 
learn python and software development. Without them, I'd hardly have found so much fun coding and
this project would never came to be. I also thank `Marcelo Reis`_ for the help naming this library :D

.. _Marcelo Reis: https://github.com/msreis
.. _Lucas Silva: https://github.com/LucasSilvaFerreira
