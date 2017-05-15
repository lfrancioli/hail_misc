
from variantqc import *
from resources import *
import re
import argparse
from rf import *
from hail.representation import *

def main(args):

    hc = hail.HailContext(log='/variantqc.log')

    if args.debug:
        logger.setLevel(logging.DEBUG)

    raw_hardcalls_path = (args.hardcalls_path if args.hardcalls_path else args.output) + ".raw_hardcalls.vds"
    raw_hardcalls_split_path = (args.hardcalls_path if args.hardcalls_path else args.output) + ".raw_hardcalls.split.vds"
    mendel_path = args.mendel_path if args.mendel_path else args.output
    rf_ann_path = args.rf_ann_path if args.rf_ann_path else args.output + ".annotated_for_rf.vds"
    rf_model_path = (args.rf_path if args.rf_path else args.output) + '.rf.model'
    rf_train_path = (args.rf_path if args.rf_path else args.output) + '.rf.training_sites.bgz'
    rf_path = args.rf_path if args.rf_path else args.output + '.rf.vds'
    rf_model = None
    vds = None


    #Create hardcalls file with raw annotations
    if args.write_hardcalls:

        variant_annotations = get_variant_type_expr()

        allele_annotations = get_stats_expr("va.stats.all_samples_raw", medians=True)
        allele_annotations.extend(get_stats_expr("va.stats.qc_samples_raw", medians=True, samples_filter_expr='sa.meta.qc_sample'))
        allele_annotations.extend(get_stats_expr("va.stats.release_samples_raw", medians=True, samples_filter_expr='sa.meta.keep'))
        allele_annotations.append("va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.nNonRefAlleles).sum()")

        (
            hc.read(full_genome_vds)
            .annotate_samples_fam(genomes_fam)
            .annotate_samples_table(hc.import_table(genomes_meta, impute=True).key_by('Sample'), root='sa.meta')
            .annotate_variants_expr(variant_annotations)
            .annotate_variants_expr("va.calldata.raw = gs.callStats(g => v) ")
            .annotate_variants_expr("va.calldata.qc_samples_raw = gs.filter(g => sa.meta.qc_sample).callStats(g => v) ")
            .annotate_variants_expr("va.calldata.release_samples_raw = gs.filter(g => sa.meta.keep).callStats(g => v) ")
            .annotate_alleles_expr(allele_annotations)
            .annotate_variants_expr(['va.nAltAlleles = v.altAlleles.filter(a => !a.isStar).length'])
            .hardcalls()
            .write(args.output + '.tmp.raw_hardcalls.vds', overwrite=True)
        )

        (
             hc.read(args.output + '.tmp.raw_hardcalls.vds')
             .min_rep()
             .repartition(num_partitions=4000)
             .write(raw_hardcalls_path, overwrite=args.overwrite)
         )

        hapmap = hc.read(hapmap_path)
        mills = hc.read(mills_path)
        omni = hc.read(omni_path)

        (
            hc.read(raw_hardcalls_path)
                .annotate_variants_expr(['va.nonsplit_alleles = v.altAlleles.map(a => a.alt)',
                                         'va.hasStar = v.altAlleles.exists(a => a.isStar)'])
                .split_multi()
                .annotate_variants_expr([
                'va.wasMixed = va.variantType == "mixed"',
                'va.alleleType = if(v.altAllele.isSNP) "snv"'
                '   else if(v.altAllele.isInsertion) "ins"'
                '   else if(v.altAllele.isDeletion) "del"'
                '   else "complex"'])
                .annotate_variants_vds(hapmap, expr='va.hapmap = isDefined(vds)')
                .annotate_variants_vds(omni, expr='va.omni = isDefined(vds)')
                .annotate_variants_vds(mills, expr='va.mills = isDefined(vds)')
                .tdt(fam=genomes_fam)
                .write(raw_hardcalls_split_path, overwrite=args.overwrite)
        )

    if args.compute_mendel:
        (
            hc.read(raw_hardcalls_split_path)
            .mendel_errors(args.output,fam=genomes_fam)
        )
        mendel_path = args.output

    #Random forests
    if args.annotate_for_rf:
        rf_ann_path = args.output + ".annotat   ed_for_rf.vds"

        rf = (
            hc.read(raw_hardcalls_split_path, drop_samples=True)
                .filter_variants_expr('va.calldata.qc_samples_raw.AC[va.aIndex] > 0')
                .annotate_variants_expr([
                'va.stats.qc_samples_raw.nrq_median = va.stats.qc_samples_raw.nrq_median[va.aIndex - 1]',
                'va.stats.qc_samples_raw.ab_median = va.stats.qc_samples_raw.ab_median[va.aIndex - 1]',
                'va.stats.qc_samples_raw.dp_median = va.stats.qc_samples_raw.dp_median[va.aIndex - 1]',
                'va.stats.qc_samples_raw.gq_median = va.stats.qc_samples_raw.gq_median[va.aIndex - 1]',
                'va.stats.qc_samples_raw.best_ab = va.stats.qc_samples_raw.best_ab[va.aIndex - 1]',
                'va.stats.qc_samples_raw.nrq = va.stats.qc_samples_raw.nrq[va.aIndex - 1]',
                'va.stats.qc_samples_raw.ab = va.stats.qc_samples_raw.ab[va.aIndex - 1]',
                'va.stats.qc_samples_raw.dp = va.stats.qc_samples_raw.dp[va.aIndex - 1]',
                'va.stats.qc_samples_raw.gq = va.stats.qc_samples_raw.gq[va.aIndex - 1]',
            ])
        )
        rf = rf.annotate_variants_table(hc.import_table(mendel_path + ".lmendel", impute=True).key_by('SNP'),
                                        expr='va.mendel = table.N')

        rf = annotate_for_random_forests(rf, hc.read(omni_path), hc.read(mills_path))
        print(rf.variant_schema)

        rf.write(rf_ann_path, overwrite=args.overwrite)

    if args.train_rf:

        features = vqsr_features if args.vqsr_features else rf_features

        if args.vqsr_training:
            tp_criteria = "va.info.POSITIVE_TRAIN_SITE"
            fp_criteria = "va.info.NEGATIVE_TRAIN_SITE"
        else:
            tp_criteria = "va.omni || va.mills"
            fp_criteria = "va.failing_hard_filters"
            if not args.no_transmitted_singletons:
                tp_criteria += " || va.transmitted_singleton"


        vds = sample_RF_training_examples(hc.read(rf_ann_path, drop_samples=True),
                                          fp_to_tp = args.fp_to_tp, tp_criteria=tp_criteria, fp_criteria=fp_criteria)

        if args.train_on_vsqsr_sites:
            vds = vds.annotate_variants_expr(['va.train = isDefined(va.info.VQSR_NEGATIVE_TRAIN_SITE) || isDefined(va.info.VQSR_POSITIVE_TRAIN_SITE)',
                                              'va.label = ...'])
        else:
            vds = vds.annotate_variants_expr(['va.train = va.train && '
                                              'va.calldata.qc_samples_raw.AC[va.aIndex] > 0 && '
                                              'v.contig != "22"'])

        rf_model = train_rf(vds, rf_features = features, num_trees =args.num_trees, max_depth = args.max_depth)
        save_model(rf_model, rf_model_path, overwrite=args.overwrite)
        vds.export_variants(rf_train_path, 'v, va.train')

    if args.apply_rf:

        if not rf_model:
            rf_model = load_model(rf_model_path)

        if not vds:
            vds = hc.read(rf_ann_path, drop_samples=True)
            vds = vds.annotate_variants_table(
                hc.import_table(rf_train_path,
                                noheader=True,
                                types = {'f0': TVariant(), 'f1':TBoolean()})
                .key_by('f0'),
                expr='va.train = table.f1')
            vds = vds.annotate_variants_expr('va.label = NA:String')

        features = vqsr_features if args.vqsr_features else rf_features

        vds = apply_rf_model(vds, rf_model, features)
        vds.write(rf_path, overwrite=args.overwrite)

    if args.write:
        #Output metrics for RF evaluation
        out_metrics = [
            'chrom = v.contig',
            'pos = v.start',
            'ref = v.ref',
            'alt = v.alt',
            'multi = va.wasSplit',
            'vqslod = va.info.VQSLOD',
            'pass = va.filters.contains("PASS")',
            'ti = v.altAllele.isTransition.toInt',
            'tv = v.altAllele.isTransversion.toInt',
            'ins = v.altAllele.isInsertion.toInt',
            'trans = va.tdt.nTransmitted',
            'untrans = va.tdt.nUntransmitted',
            'type = va.variantType',
            'qd = va.info.QD',
            'train_vqsr = va.info.NEGATIVE_TRAIN_SITE || va.info.POSITIVE_TRAIN_SITE',
            'label_vqsr = if(va.info.POSITIVE_TRAIN_SITE) "TP" else orMissing(isDefined(va.info.NEGATIVE_TRAIN_SITE), "FP")',
            'ac_origin = va.info.AC[va.aIndex-1]',
            'an_origin = va.info.AN',
            'ac_unrelated = va.AC_unrelated[va.aIndex-1]',
            'mendel_err = va.mendel'
        ]

        rf_out = (
            hc.read(rf_ann_path, drop_samples=True)
            .filter_variants_table(hail.KeyTable.import_interval_list(lcr_path), keep=False)
            .filter_variants_table(hail.KeyTable.import_interval_list(decoy_path), keep=False)
            .filter_variants_expr('va.calldata.qc_samples_raw.AC[va.aIndex] > 0')
        )

        if args.rf_ann_files:
            rf_out, additional_metrics = annotate_with_additional_rf_files(rf_out,args.rf_ann_files)
            out_metrics.extend(additional_metrics)

        #Get number of SNVs and indels to sample
        nvariants = rf_out.query_variants(["variants.filter(x => x.altAllele.isSNP).count()",
                               "variants.filter(x => x.altAllele.isIndel).count()"])

        logger.info("Number of SNVs: %d, number of Indels: %d" % (nvariants[0], nvariants[1]))

        (
            rf_out
            .annotate_variants_table(hc.import_table(mendel_path + ".lmendel", impute=True).key_by('SNP'),expr='va.mendel = table.N')
            .filter_variants_expr('(v.altAllele.isSNP && pcoin(2500000.0 / %d) )|| '
                                  '(v.altAllele.isIndel && pcoin(2500000.0 / %d) )||'
                                  '(va.mendel>0 && va.calldata.raw.AC[va.aIndex] == 1 )' % (nvariants[0],nvariants[1]))
            .export_variants(rf_path + ".va.txt.bgz", ",".join(out_metrics))
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write_hardcalls', help='Creates a hardcalls vds', action='store_true')
    parser.add_argument('--hardcalls_path', help='Overrides the default hardcalls paths prefix ($output)')
    parser.add_argument('--compute_mendel', help='Computes Mendel errors', action='store_true')
    parser.add_argument('--mendel_path', help='Overrides the default mendel path prefix ($output)')
    parser.add_argument('--annotate_for_rf', help='Creates an annotated VDS with features for RF', action='store_true')
    parser.add_argument('--rf_ann_path', help='Overrides the default rf annotation path ($output + .annotated_for_rf.vds)')
    parser.add_argument('--train_rf', help='Trains RF model', action='store_true')
    parser.add_argument('--rf_model_path',
                        help='Overrides the default rf model prefix ($output)')
    parser.add_argument('--rf_path',
                        help='Overrides the default rf path ($output + .rf.vds)')
    parser.add_argument('--apply_rf', help='Applies RF model to the data', action='store_true')
    parser.add_argument('--write', help='Writes the results of the current RF model with the given name.', action='store_true')
    parser.add_argument('--rf_ann_files', help='RF files to annotate results with in pipe-delimited format: name|location|rf_root', nargs='+')
    parser.add_argument('--fp_to_tp', help='Sets to ratio of TPs to FPs for creating the RF model.', default=1.0, type=float)
    parser.add_argument('--num_trees', help='Number of tress in the RF model.', default=500)
    parser.add_argument('--max_depth', help='Maxmimum tree depth in the RF model.', default=5)
    parser.add_argument('--vqsr_features', help='Use VQSR features only (+ snv/indel variant type)', action='store_true')
    parser.add_argument('--vqsr_training', help='Use VQSR training examples',
                        action='store_true')
    parser.add_argument('--no_transmitted_singletons', help='Do not use transmitted singletons for training.',
                        action='store_true')
    parser.add_argument('--train_on_vsqsr_sites', help='Use VQSR training sites',
                        action='store_true')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    args = parser.parse_args()
    main(args)
