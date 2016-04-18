import numpy as np


class UnsupportedFile(Exception):
    pass


def detect_mime(path_to_file, uncompress=False):
    import magic
    mime = magic.Magic(mime=True, uncompress=uncompress)
    mime_from_file = mime.from_file(path_to_file)
    # magic returns byte, not str
    # ToDo: remove decode when issue #98 is closed
    #  https://github.com/ahupp/python-magic/issues/98
    return mime_from_file.decode()


def magic_open(path_to_file, mode='rt'):
    import gzip
    mime = detect_mime(path_to_file)

    if mime == 'text/plain':
        return open(path_to_file, mode=mode)

    elif mime == 'application/x-gzip':
        # if detect_mime(path_to_file, uncompress=True) == 'text/plain':
        return gzip.open(path_to_file, mode=mode)
    else:
        raise UnsupportedFile('File %s is type %s' % (path_to_file, mime))


array2str = lambda arr: ','.join(str(x) for x in arr)
str2array = lambda string: np.fromstring(string, sep=',', dtype=int)


def rand_id():
    import string as s
    chars = s.ascii_lowercase + s.digits
    chars = [x for x in chars]
    _id = ''.join(np.random.choice(chars, 15))
    return "Random_ID_" + _id
