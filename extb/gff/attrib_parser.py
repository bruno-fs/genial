import re

from ..exceptions import ParseError, UnsupportedFile


def attributes_parser(attributes: str, file_format='gff3') -> dict:
    def compile_pattern(ff):
        if ff == 'gff3':
            return re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
        elif ff == 'gtf':
            return re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
        else:
            raise UnsupportedFile('Unsupported file_format: %s' % ff)

    # Workarounds for incorrect use of SEMICOLONS
    # ToDo: implement a function for this? is it better to check for each case?
    """
    SEMICOLON is a SACRED character on the attributes field of GFF/GTF files
            they should ONLY be used to separate attributes

    These workarounds will remove known cases of misuses of semicolons
    """

    # 1: escaped html characters
    """
    Some gffs from WormBase have escaped HTML characters like &amp; (&)
    unescaping will get us rid of them

    However, there is a html code for SEMICOLON and we don't want to unescape it.
    Let's replace it by its percent encoded version ;)

    character |  html escaped | percent escaped
    --------- | ------------- | ---------------
        ;     |     &#59;     |      %3B


    """
    from html import unescape
    attributes.replace("&#59;", "%3B")
    attributes = unescape(attributes)

    # 2: semicolons inside atributes! (What's wrong with u, ensembl???)
    """
    Ensembl 2014 Homo_sapiens.GRCh38.78.gtf has* lines containing attributes like:

    gene_name "PRAMEF6;"
    transcript_name "PRAMEF6;-201"

    This is SO UGLY and WRONG. Let's simple remove those atrocities.

    Of course we could percent-escape'em, but if these semicolons were so important to be included,
    the author of the GFF would have escaped them in first place.

    source: https://gist.github.com/iansealy/b9cbc56bd1affe10d37a#count-reads-with-htseq-for-tophat-bam-files
    * it seems ensembl fixed and replaced those files on ftp in 2015, but some people might
    still be using the f*ck*d up version provided in 2014.

    """

    attributes.replace(';\"', '\"').replace(";-", "-")

    attrib_dict = {}
    attribs = re.sub(';\s*$', '', attributes)
    attribs = attribs.split(';')
    pattern = compile_pattern(file_format)

    for att in attribs:
        g = re.search(pattern, att)

        from sys import intern
        # intern strings to save memory
        try:
            k, v = g.group(1, 2)

            if file_format == 'gff3' and re.match(r'ID|Parent', k):
                has_ids_prepended = re.match(r'^(transcript|gene):', v)
                if has_ids_prepended:
                    # ensembl GFF usually has transcript/gene: prepended
                    # to the values of ID/Parent attributes
                    # replace transcript/gene: prepended to values
                    v = re.sub(r'^(transcript|gene):', '', v)

                    id_prepended = has_ids_prepended.group(1) + '_id'

                    # let's use this prepended information in our advantage
                    # and create an gene_id/transcript_id for GFF files!

                    if id_prepended not in attrib_dict:
                        attrib_dict[intern(id_prepended)] = intern(v)

            attrib_dict[intern(k)] = intern(v)
        except AttributeError:
            # TODO: use/make more specific exceptions
            raise ParseError('regex for %s failed to parse %s' % (file_format, attributes))

    return attrib_dict