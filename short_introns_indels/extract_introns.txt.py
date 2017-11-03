import argparse
import gzip

def print_introns(input, output):
    for line in input:
        if not line.startswith("#"):
            fields = line.split("\t")
            exons = zip(fields[9].rstrip(",").split(","),
                        fields[10].rstrip(",").split(",")
                        )

            if len(exons) != int(fields[8]):
                print "WARN: {} exons found for transcript {}, but exonCount col = {}".format(
                    len(exons),
                    fields[1],
                    fields[8]
                )

            last_end = 0
            for start, end in exons:
                start = int(start)
                end = int(end)
                if last_end > 0 and last_end+1 < start-1:
                    output.write("\t".join([
                        fields[2][3:] if fields[2].startswith("chr") else fields[2],
                        str(last_end + 1),
                        str(start - 1),
                        fields[1]

                    ]) + "\n"
                )
                last_end = end


def main(args):

    with open(args.output, 'w') as output:
        if args.input.endswith(".gz"):
            with gzip.open(args.input, 'rb') as input:
                print_introns(input, output)
        else:
            with open(args.input) as input:
                print_introns(input, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Input UCSC file', required=True)
    parser.add_argument('--output', help='Output BED file', required=True)
    args = parser.parse_args()

    main(args)
