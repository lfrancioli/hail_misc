from hail import *
from resources import *
import argparse
from utils import *

def main(args):
    hc = HailContext()

    #logger.info("Loading introns from bed and creating interval list")

    ends_length = args.ends_length

    introns_kt = hc.import_table(args.bed, no_header=True, impute=True, types={'f0': TString()})
    introns_kt = introns_kt.repartition(20)
    introns_kt = introns_kt.rename({'f0': 'chrom', 'f1': 'start', 'f2': 'end', 'f3': 'transcript', 'f4': 'sequence'})
    introns_kt = introns_kt.annotate(['interval = Interval(Locus(chrom, start), Locus(chrom, end))',
                                      'metric_entropy = sequence.entropy() / sequence.length',
                                      'start_entropy = sequence[:{0}].entropy() / sequence[:{0}].length'.format(ends_length),
                                      'end_entropy = sequence[-{0}:].entropy() / sequence[-{0}:].length'.format(ends_length),
                                      'ends_entropy = let s = sequence[-{0}:] + sequence[:{0}] in s.entropy() / s.length'.format(ends_length),
                                      'middle_entropy = sequence[{0}:-{0}].entropy() / sequence[{0}:-{0}].length'.format(ends_length),
                                      'intron_id = transcript + "_" + str(start)'])
    introns_kt = introns_kt.persist()

    introns_intervals = introns_kt.query('interval.collect()')

    vds = get_gnomad_public_data(hc, "genomes", split=True)
    vds = vds.filter_intervals(introns_intervals)
    vds = vds.filter_variants_expr("v.altAllele().isIndel()")
    vds = vds.annotate_variants_table(introns_kt.select(['interval','intron_id','start','end']).key_by(['interval']),
                                      root = 'va.intron')

    indel_intron_counts_kt = vds.variants_table()

    indel_intron_counts_kt = indel_intron_counts_kt.annotate([
        'distance_to_start = v.start - va.intron.start',
        'distance_to_end = va.intron.end - v.start'
    ])
    indel_intron_counts_kt = indel_intron_counts_kt.aggregate_by_key(['intron_id = va.intron.intron_id',
                                                                      'deletion = v.altAllele.isDeletion()',
                                                                      'pass = va.filters.isEmpty()',
                                                                      'lcr = va.info.lcr',
                                                                      'segdup = va.info.segdup',
                                                                      'wasSplit = va.wasSplit'],
                                                                     [
                                                                         'n = v.count()',
                                                                         'n_start = distance_to_start.filter(x => x <= {}).count()'.format(ends_length),
                                                                         'n_end = distance_to_end.filter(x => x <= {0}).count()'.format(ends_length)
                                                                     ])
    indel_intron_counts_kt = indel_intron_counts_kt.key_by('intron_id')


    introns_kt = introns_kt.key_by(['intron_id']).join(indel_intron_counts_kt, how="left")
    introns_kt.drop(['interval','sequence']).export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed', help='Input introns BED file', required=True)
    parser.add_argument('--out', help='Output file', required=True)
    parser.add_argument('--ends_length', help='Length of the window around ends to aggregate', default = 100, required=False, type=int)
    args = parser.parse_args()

    main(args)