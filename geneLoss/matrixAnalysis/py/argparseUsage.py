# how to use argparse?
# https://docs.python.org/3/library/argparse.html


import argparse

# I - The first step in using the argparse is creating an ArgumentParser object:
parser = argparse.ArgumentParser(description='Process some integers.')

# II - Filling an ArgumentParser with information about program arguments is done by making calls to the add_argument() method.
parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')


