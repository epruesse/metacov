import sys

if sys.version_info[0] > 2:
    FileNotFoundError = FileNotFoundError
else:
    FileNotFoundError = IOError
