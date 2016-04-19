class ParseError(Exception):
    pass


class UnsupportedFile(Exception):
    pass


class MultipleParentsGFF(UnsupportedFile):
    pass

