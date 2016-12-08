import re
from warnings import warn
from genial.exceptions import ParseError, UnsupportedFile


def attributes_parser(attributes: str, file_format='gff3') -> dict:
    """
     SEMICOLON is a SACRED character on the attributes field of GFF/GTF files
             they should ONLY be used to separate attributes

     These workarounds will remove known cases of misuses of semicolons

    # 1: escaped html characters

    Some gffs from WormBase have escaped HTML characters like &amp; (&)
    unescaping will get us rid of them

    However, there is a html code for SEMICOLON and we don't want to unescape it.
    Let's replace it by its percent encoded version ;)

    character |  html escaped | percent escaped
    --------- | ------------- | ---------------
        ;     |     &#59;     |      %3B


    # 2: semicolons inside atributes! (What's wrong with u, ensembl???)

    Ensembl 2014 Homo_sapiens.GRCh38.78.gtf has* lines containing attributes like:

    gene_name "PRAMEF6;"
    transcript_name "PRAMEF6;-201"

    This is SO UGLY and WRONG. Let's simple remove those atrocities.

    Of course we could percent-escape'em, but if these semicolons were so important to be included,
    the author of the GFF would have escaped them in first place.

    source: https://gist.github.com/iansealy/b9cbc56bd1affe10d37a#count-reads-with-htseq-for-tophat-bam-files
    * it seems ensembl fixed and replaced those files on ftp in 2015, but some people might
    still be using the f* up version provided in 2014.


     """
    def compile_pattern(ff):
        if ff == 'gff3':
            return re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
        elif ff == 'gtf':
            return re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
        # else:
        #     raise UnsupportedFile('Unsupported file_format: %s' % ff)

    # Workarounds for incorrect use of SEMICOLONS
    # ToDo: implement a function for this? is it better to check for each case?
    # 1: escaped html characters
    from html import unescape
    attributes.replace("&#59;", "%3B")
    attributes = unescape(attributes)

    # 2: semicolons inside atributes! (What's wrong with u, ensembl???)
    attributes.replace(';\"', '\"').replace(";-", "-")

    # 3: trailing semicolon
    attributes = re.sub(';\s*$', '', attributes)

    # finally, get a list of attributes
    attrib_list = attributes.split(';')

    from genial.utils import InternDict
    attrib_dict = InternDict()
    # attrib_dict = {}

    # compile the pattern to parse_to_dict gff3 or gtf
    pattern = compile_pattern(file_format)
    for att in attrib_list:
        g = re.search(pattern, att)

        try:
            k, v = g.group(1, 2)
        except AttributeError:
            # TODO: use/make more specific exceptions
            # raise ParseError('regex for %s failed to parse_to_dict %s' % (file_format, attributes))
            warn('regex for %s failed to parse attribute %s' % (file_format, att))
        # if file_format == 'gff3' and re.match(r'ID|Parent', k):
        #     has_ids_prepended = re.match(r'^(transcript|gene):', v)
        #     if has_ids_prepended:
        #         # ensembl GFF usually has transcript/gene: prepended
        #         # to the values of ID/Parent attributes
        #         # replace transcript/gene: prepended to values
        #         v = re.sub(r'^(transcript|gene):', '', v)
        #
        #         id_prepended = has_ids_prepended.group(1) + '_id'
        #
        #         # let's use this prepended information in our advantage
        #         # and create an gene_id/transcript_id for GFF files!
        #
        #         if id_prepended not in attrib_dict:
        #             attrib_dict[id_prepended] = v
        else:
            has_ids_prepended = re.match(r'^(\w+):', v)
            if has_ids_prepended:
                v = re.sub(r'^\w+:', '', v)

            attrib_dict[k] = v

    return attrib_dict