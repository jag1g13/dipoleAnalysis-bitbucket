import sys


class FrameReader(object):
    def __init__(self):
        print("Constructed FrameWriter")

    def __enter__(self, filename):
        self.file = open(filename, 'w')
        self.filename = filename
        print("Opened file: {0}".format(self.filename))
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()
        print("Closed file: {0}".format(self.filename))

if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)
