import hail
import argparse
from utils import *


def main(args):

    hc = hail.HailContext(log = "/sampleqc.log")

    if args.import_vcf:
        vds = hc.import_vcf(args.import_vcf)
        vds = vds.sampleqc()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--import_vcf', help='VCF file to import. Imports the given vcf file and annotates it with VEP and sampleQC metrics')
    parser.add_argument('--export_sample_qc', help='VDS file to get sampleQC metrics from. Exports the sampleQC metrics')
    parser.add_argument('--run_mendel', help='VDS file to run mendel_errors on. Runs the sampleQC on the given file.')
    parser.add_argument('--output', help='Output path.')



