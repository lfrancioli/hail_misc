from variantqc import *
import re

try:
    hc
except NameError:
    hc = hail.HailContext(log='/variantqc.log')


# Inputs
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats22.vds"
iteration_re = r"newStats(\d+)"
other_rf_ann_files = ["gs://gnomad/RF/gnomad.sites.RF.newStats17.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats18.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats19.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats20.vds"]
raw_hardcalls_split_path = "gs://gnomad/gnomad.raw_hardcalls.split.vds"
exomes_hardcalls_path = "gs://exac2/hardcalls/exacv2.raw.hardcalls.splitmulti.qc.vds"
exomes_to_combined_IDs = "gs://gnomad/exac_to_combined.IDs.txt"
exomes_concordance_samples = "gs://gnomad/exac_samples_in_gnomad.txt"
genomes_to_combined_IDs = "gs://gnomad/gnomad_to_combined.IDs.txt"
genomes_concordance_samples = "gs://gnomad/gnomad_samples_in_exac.txt"
# Missing
# exomes_high_conf_regions_path = "gs://exac2/high_coverage.auto.interval_list"

# Outputs
syndip_concordance_prefix = "gs://gnomad/truth-sets/gnomad_hybrid"
NA12878_concordance_prefix = "gs://gnomad/truth-sets/gnomad_NA12878"
exomes_concordance_prefix = "gs://gnomad/gnomad_exac_raw"

# Annotations to output
concordance_annotations = ['chrom = v.contig',
'pos = v.start',
'ref = v.ref',
'alt = v.alt',
'snp = v.altAllele.isSNP',
'ins = v.altAllele.isInsertion',
'del = v.altAllele.isDeletion',
'concordance = va.concordance'
                           ]

rf_features_annotations = [
                'va.rf.info.QD',
                'va.rf.info.MQ',
                'va.rf.info.MQRankSum',
                'va.rf.info.FS',
                'va.rf.info.SOR',
                'va.rf.info.InbreedingCoeff',
                'va.rf.info.ReadPosRankSum',
                'va.rf.stats.raw.nrq_median',
                'va.rf.stats.raw.ab_median',
                'va.rf.stats.raw.dp_median',
                'va.rf.stats.raw.gq_median']

truth_concordance_annotations = list(concordance_annotations)
out_num = re.search(iteration_re,rf_path).groups()[0]
truth_concordance_annotations.extend(['type = va.variantType',
                                      'wassplit = va.left.wasSplit',
                                      'vqslod = va.rf.info.VQSLOD',
                                      'truth.wassplit = va.right.wasSplit',
                                      'truth_gt = va.truth_gt',
                                      'called_gt = va.called_gt',
                                      'training = va.rf.train',
                                      'label = va.rf.label',
                                      'rfpred%s = va.rf.RF1.prediction' % out_num,
                                      'rfprob%s = va.rf.RF1.probability["TP"]' % out_num
                                      ])

exomes_concordance_annotations = list(concordance_annotations)
exomes_concordance_annotations.extend(['gnomad.multi = va.left.wasSplit',
                                       'gnomad.vqslod = va.left.info.VQSLOD',
                                       'gnomad.qd = va.left.info.QD',
                                       'gnomad.ac = va.left.info.AC[va.right.aIndex-1]',
                                       'gnomad.an = va.left.info.AN',
                                       'exac.multi = va.right.wasSplit',
                                       'exac.vqslod = va.right.info.VQSLOD',
                                       'exac.ac = va.right.info.AC[va.left.aIndex-1]',
                                       'exac.an = va.right.info.AN',
                                       'exac.qd = va.right.info.QD',
                                       'gnomad.training = va.rf.train',
                                       'gnomad.label = va.rf.label',
                                       'gnomad.rfpred1 = va.rf.RF1.prediction',
                                       'gnomad.rfprob1 = va.rf.RF1.probability["TP"]'
                                       ])

for f in other_rf_ann_files:
    out_num = f[39:41]
    truth_concordance_annotations.append(
        'rfpred%s = va.rf.RF%s.prediction,rfprob%s = va.rf.RF%s.probability["TP"]' % (out_num, out_num, out_num, out_num))

print(truth_concordance_annotations)

# Actions
compute_syndip_concordance = False
export_syndip_concordance = True
compute_NA12878_concordance = False
export_NA12878_concordance = True
compute_Exomes_concordance = False
export_Exomes_concordance = False

def add_annotations_for_concordance(vds):
    for f in other_rf_ann_files:
        out_num = re.search(iteration_re,f).groups()[0]
        vds = vds.annotate_variants_vds(hc.read(f, sites_only=True), code=
                                              'va.RF%s.prediction = vds.RF1.prediction,'
                                              'va.RF%s.probability = vds.RF1.probability' % (out_num,out_num))
    return vds

if compute_syndip_concordance:
    compute_concordance(hc.read(raw_hardcalls_split_path),
                        hc.read(syndip_path),
                        'CHMI_CHMI3_WGS1',
                        syndip_high_conf_regions_path,
                        syndip_concordance_prefix)

if export_syndip_concordance:
    rf_vds = add_annotations_for_concordance(hc.read(rf_path, sites_only=True).min_rep())
    add_annotations_for_concordance(rf_vds)
    export_concordance(hc.read(syndip_concordance_prefix + ".v_concordance.vds"),
                       rf_vds,
                       truth_concordance_annotations,
                       syndip_concordance_prefix)


if compute_NA12878_concordance:
    compute_concordance(hc.read(raw_hardcalls_split_path),
                        hc.read(NA12878_path).rename_samples("gs://gnomad/truth-sets/NA12878.gnomad_rename"),
                        'G94982_NA12878',
                        NA12878_high_conf_regions_path,
                        NA12878_concordance_prefix)

if export_NA12878_concordance:
    rf_vds = add_annotations_for_concordance(hc.read(rf_path, sites_only=True).min_rep())
    export_concordance(hc.read(NA12878_concordance_prefix + ".v_concordance.vds"),
                       rf_vds,
                       truth_concordance_annotations,
                       NA12878_concordance_prefix)


if (compute_Exomes_concordance or export_Exomes_concordance):
    if(compute_Exomes_concordance):
        exomes = (
            filter_for_concordance(hc.read(exomes_hardcalls_path), high_conf_regions=exomes_high_conf_regions_path)
            .filter_samples_list(exomes_concordance_samples)
            .rename_samples(exomes_to_combined_IDs)
            .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count()>0', keep=True)
        )
        exomes.annotate_global_expr_by_sample('global.nsamples = samples.count()').show_globals()

        genomes = (
            filter_for_concordance(hc.read(raw_hardcalls_split_path), high_conf_regions=exomes_high_conf_regions_path)
            .filter_samples_list(genomes_concordance_samples)
            .rename_samples(genomes_to_combined_IDs)
            .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count()>0', keep=True)
        )
        genomes.annotate_global_expr_by_sample('global.nsamples = samples.count()').show_globals()

        (s_concordance, v_concordance) = genomes.concordance(right=exomes)

        s_concordance.write(exomes_concordance_prefix + ".s_concordance.vds")
        v_concordance.write(exomes_concordance_prefix + ".v_concordance.vds")

    else:
        v_concordance = hc.read(exomes_concordance_prefix + ".v_concordance.vds")

    if(export_Exomes_concordance):
        (
            v_concordance.annotate_variants_vds(hc.read(rf_path), root='va.rf')
                .filter_variants_expr('va.concordance[3].exists(x => x>0) || va.concordance[4].exists(x => x>0) || va.concordance[0:2].map(x => x[3:4].sum).sum >0')
            .export_variants(exomes_concordance_prefix + ".stats.txt.bgz", ",".join(exomes_concordance_annotations))
        )