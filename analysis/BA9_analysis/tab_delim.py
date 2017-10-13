#!/usr/bin/env python
import pandas as pd
import sys

def readFile(filename):
    """Read file"""
    df = pd.read_csv(filename, sep="\t")
    return df


def main():
    if len(sys.argv) < 1:
        print("Usage: {0} <CSV File to be tab-delimited>".format(sys.arv[0]))
        sys.exit(1)

    file1 = sys.argv[1]
    df = readFile(file1)
    df.to_csv("BA9_salmon_firth.csv", index=False)

if __name__ == "__main__":
    main()
