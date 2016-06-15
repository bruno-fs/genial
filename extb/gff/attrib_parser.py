import re

from extb.exceptions import ParseError


def attributes_parser(attributes, file_format='gff3'):
    def compile_pattern(ff):
        if ff == 'gff3':
            return re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
        elif ff == 'gtf':
            return re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
        else:
            raise Exception('Unsupported file_format: %s' % ff)

    from html import unescape
    # unescape HTML characters like &amp;
    # necessary on some gffs from WormBase
    attributes = unescape(attributes)

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

            if re.match(r'ID|Parent', k):
                has_ids_prepended = re.match(r'^(transcript|gene):', v)
                if has_ids_prepended:
                    # ensembl GFF usually has transcript/gene: prepended
                    # to the values of ID/Parent attributes
                    # replace transcript/gene: prepended to values
                    # v = re.sub(r'^(transcript|gene):', '', v)
                    v = re.sub(r'^(\w+):', '', v)

                    id_prepended = has_ids_prepended.group(1) + '_id'
                    import sys
                    # print(id_prepended, v, file=sys.stderr)
                    # let's use this prepended information in our advantage
                    # and create an gene_id or transcript_id for the attributes
                    if id_prepended not in attrib_dict:
                        attrib_dict[intern(id_prepended)] = intern(v)

            attrib_dict[intern(k)] = intern(v)
        except AttributeError:
            # TODO: use/make more specific exceptions
            raise ParseError('regex for %s failed to parse %s' % (file_format, attributes))

    return attrib_dict