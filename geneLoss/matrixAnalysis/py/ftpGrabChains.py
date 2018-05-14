#!/usr/bin/env python

import subprocess
import sys

USAGE = '''Usage: %s hg19 assemblyList

Will download chains from target to all assemblies in assemblyList
Assumes assemblies are lowercase (e.g. hg19, bostau7)

Note: will fail on some machines (e.g. ttn) because lftp is not installed.
Runs on dev
''' % sys.argv[0]

HOST='hgdownload.cse.ucsc.edu'
USER='anonymous'
PASS='robinjia@stanford.edu'

def makeGetString(target, query):
    queryUpper = query[0].upper() + query[1:]
    return ("get goldenPath/%s/vs%s/%s.%s.all.chain.gz"
                % (target, queryUpper, target, query))

def main():
    target = sys.argv[1]
    asmFile = sys.argv[2]
    with open(asmFile) as f:
        getStrings = [ makeGetString(target, line.strip()) for line in f ]
    getStrings.append('exit;')
    getCommand = '\n'.join(getStrings)

    try:
        subprocess.check_call(["lftp", "-e", getCommand,
                "-u", "%s,%s" % (USER, PASS), HOST])
    except Exception as e:
        print >>sys.stderr, "Error in calling lftp (%d): %s" % (e.errno,
                e.strerror)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print >> sys.stderr, USAGE
        exit(1)
    main()
