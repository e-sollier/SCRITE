import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--minREF", type=int, default=200,help='Minimum number of reference reads.')
parser.add_argument("--minALT", type=int, default=30,help='Minimum number of alternative reads.')

args = parser.parse_args()

for line in sys.stdin:
    if line[0]=="#":
        sys.stdout.write(line)
    else:
        linesplit = line.split("\t")
        info = linesplit[9].split(":")
        n_reads_ref = int(info[4])
        n_reads_alt = int(info[5])
        if n_reads_ref>=args.minREF and n_reads_alt>=args.minALT: 
            sys.stdout.write(line)

