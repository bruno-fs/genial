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
    attributes = unescape(attributes)

    attrib_dict = {}
    attribs = re.sub(';\s*$', '', attributes)
    attribs = attribs.split(';')
    pattern = compile_pattern(file_format)

    for att in attribs:
        g = re.search(pattern, att)

        try:
            k, v = g.group(1, 2)
            # replace transcript:/gene: prepended to these values
            # like
            v = re.sub(r'^(transcript|gene):', '', v)
            attrib_dict[k] = v
        except AttributeError:
            # TODO: use/make more specific exceptions
            raise ParseError('regex %s failed to parse %s' % (str(pattern), attributes))

    return attrib_dict