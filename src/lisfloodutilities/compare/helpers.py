import argparse
import sys


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):

        self.add_argument("-a", "--dataset_a", help='path to outputh of LisFlood version A', required=True)
        self.add_argument("-b", "--dataset_b", help='path to outputh of LisFlood version B', required=True)
        self.add_argument("-m", "--maskarea", help='path to mask', required=True)
