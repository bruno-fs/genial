genial: GENome Interactive Annotation Library
=============================================

.. image:: https://badge.fury.io/py/genial.png
    :target: https://badge.fury.io/py/genial

This library provides `InteractiveAnnotation`, a high level representation of genome annotations,
allowing users to easily extract information, manipulate and reformat commonly used annotation
formats such as GFF and BED.

The package is currently in alpha stage and only runs on `python3`.

For convenience, we provided two CLI utilities:

.. code-block:: bash

    $ annotParser.py -h
    usage: annotParser.py [-h] [-i INPUT] [-o OUTPUT] [-f INPUT_FORMAT]
                          [-t OUTPUT_FORMAT] [-n MIN_EXON_COUNT]
                          [-igs IGNORE_GAPS_SMALLER_THAN]
                          [-igb IGNORE_GAPS_BIGGER_THAN] [-v]

    parse, filter and convert annotation files

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            input file. to read from pipe, use the argument
                            'stdin'
      -o OUTPUT, --output OUTPUT
                            output file
      -f INPUT_FORMAT, --input_format INPUT_FORMAT
                            input file format (supported formats: gtf, gff3, bed)
      -t OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                            output file format (supported formats: extb, bed)
      -n MIN_EXON_COUNT, --min_exon_count MIN_EXON_COUNT
                            min number of exons (inclusive)
      -igs IGNORE_GAPS_SMALLER_THAN, --ignore_gaps_smaller_than IGNORE_GAPS_SMALLER_THAN
      -igb IGNORE_GAPS_BIGGER_THAN, --ignore_gaps_bigger_than IGNORE_GAPS_BIGGER_THAN
      -v, --invert_match    only return false results (similar to grep -v)


.. code-block:: bash

    $ annotMergeSmallGaps.py -h
    usage: annotMergeSmallGaps.py [-h] [-i INPUT] [-o OUTPUT] [-f INPUT_FORMAT]
                                  [-t OUTPUT_FORMAT] [-s SMALL_GAP_SIZE]

    Merge exons separated by small gaps. Can also be used to convert different
    kinds of annotations.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            input file. to read from pipe, use the argument
                            'stdin'
      -o OUTPUT, --output OUTPUT
                            output file
      -f INPUT_FORMAT, --input_format INPUT_FORMAT
                            input file format (supported formats: gff3, bed, gtf)
      -t OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                            output file format (supported formats: extb, bed)
      -s SMALL_GAP_SIZE, --small_gap_size SMALL_GAP_SIZE
                            gap size. exons separated by gaps with this size or
                            less should be removed



Both the scripts above can also be used to convert from different kinds of annotation files.
A more advanced usage can be achieved importing the library.

Supported formats
-----------------

* Input: GFF3, GTF, bed

* Output: bed


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
