import argparse
import subprocess
import os

def main(args):
    chrom = ""
    pos = 0

    with open(args.sorted_partitions) as partitions_file:
        for p in partitions_file:
            c, s, e = p.strip().split("\t")
            if c == chrom and int(s) <= pos:
                print "{}\t{}\t{}".format(c,pos,s)

                with open(args.vcf_starts) as vcf_files:
                    for l in vcf_files:
                        vcf = l.strip().split("\t")
                        if len(vcf) > 2 and vcf[1] == c and vcf[2] == s:
                            print "{}\t{}\t{}".format(vcf[0], vcf[1], vcf[2])

                            vcf_name = os.path.split(vcf[0])[1]
                            # tabix_commands = [
                            #     "gsutil cp {} .".format(vcf[0]),
                            #     "gsutil cp {}.tbi .".format(vcf[0]),
                            #     "gsutil mv {0} {0}.old".format(vcf[0]),
                            #     "tabix -h {0} {1}:{2}-300000000 | bgzip -c > {0}.new".format(vcf_name, vcf[1], pos + 1),
                            #     "gsutil cp {}.new {}".format(vcf_name, vcf[0]),
                            #     "rm {0} {0}.tbi {0}.new".format(vcf_name)
                            # ]

                            commands = [
                                "gsutil cp {}.old ./{}".format(vcf[0], vcf_name),
                                "gsutil cp {}.tbi .".format(vcf[0]),
                                "cat <(tabix -H {0}) <(zcat {0} | awk '/^[^#]/ {{FS=\"\\t\"; if($1 != \"{1}\" || $2 > {2}){{print;}}}}') | bgzip -c > {0}.new".format(vcf_name, vcf[1], pos),
                                "gsutil cp {}.new {}".format(vcf_name, vcf[0]),
                                "rm {0} {0}.tbi {0}.new".format(vcf_name)
                            ]
                            for command in commands:
                                if args.print_only:
                                    print command
                                else:
                                    subprocess.check_output(command, shell=True, executable='/bin/bash')
            chrom = c
            pos = int(e)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_starts', help='Files containing start positions of all VCF files')
    parser.add_argument('--sorted_partitions', help='File containing sorted partitions')
    parser.add_argument('--print_only', help='Does not run anything -- just prints commands.', action='store_true')
    args = parser.parse_args()
    main(args)