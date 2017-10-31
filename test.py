
from utils import *
import hail
from hail.representation import Interval
from pprint import pprint

#Inputs
vds_path = "gs://gnomad/gnomad.10ksites.vds"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats24.vds"
vep_config = "/vep/vep-gcloud.properties"
mendel_path = "gs://gnomad/gnomad.raw_calls"
fam_path = "gs://gnomad/gnomad.final.goodTrios.fam"
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"

genome_pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']

hc = hail.HailContext(log='/test.log')
#hc.write_partitioning(full_exome_01_vds_path) ##TODO: Remove this

def get_ab_stats(outprefix):
    for interval, filter_males in [("1-22",False),("X",True)]:
        vds = get_gnomad_data(hc, "genomes", release_samples=True)
        vds = vds.filter_intervals(Interval.parse(interval))
        if filter_males:
            vds = vds.filter_samples_expr('sa.meta.sex == "female"')
        vds = vds.split_multi()
        products = vds.query_samples('samples.map(s => sa.meta.Product).collectAsSet()')

        products_dict = {product: re.sub('\s+','_',re.sub('[^0-9a-zA-Z\s]+', '', product)).lower() for product in products}

        annotations_desc ={
            'HETC': {'Description': 'Total number of heterozygous genotypes'},
            'REFC': {'Description':'Total number of reference alleles across heterozygous genotypes'},
            'ALTC': {'Description':'Total number of alternate alleles across heterozygous genotypes'},
            'ALTR': {'Description':'Ratio of alternate alleles across heterozygous genotypes'},
            'ARP': {'Description':'P-value for imbalance of alternate allele ratio across heterozygous genotypes'}
        }

        expr = []
        for p, pname in products_dict.iteritems():
            expr.extend([
                'va.info.{0}.HETC = gs.filter(g => sa.meta.Product == "{1}" && g.isHet).count()'.format(pname, p),
                'va.info.{0}.REFC = gs.filter(g => sa.meta.Product == "{1}" &&  g.isHet).map(g => g.ad[0]).sum()'.format(pname, p),
                'va.info.{0}.ALTC = gs.filter(g => sa.meta.Product == "{1}" &&  g.isHet).map(g => g.ad[1]).sum()'.format(pname, p),
                'va.info.{0}.ALTR = let all_ad = gs.filter(g => sa.meta.Product == "{1}" &&  g.isHet).map(g => g.ad).sum() in all_ad[1]/all_ad.sum()'.format(
                    pname, p),
                'va.info.{0}.ARP = let all_ad = gs.filter(g => sa.meta.Product == "{1}" &&  g.isHet).map(g => g.ad).sum() in binomTest(all_ad[1], all_ad.sum(), 0.5, "two.sided")'.format(
                    pname, p)
            ])

        vds = vds.annotate_variants_expr(expr)
        vds = vds.drop_samples()
        #vds = vds.repartition(50 if filter_males else 1000)
        vds.write("{}.chr{}.vds".format(outprefix, interval))

        vds = hc.read("{}.chr{}.vds".format(outprefix, interval))
        for pname in products_dict.values():
            pvds = vds.annotate_variants_expr('va.info = va.info.{}'.format(pname))
            pvds = pvds.filter_variants_expr('va.info.HETC > 0')
            pvds = pvds.annotate_variants_expr(['va.info = select(va.info, {})'.format(",".join(annotations_desc.keys())),
                                                'va.filters = NA:Set[String]',
                                                'va.rsid = NA:String',
                                                'va.qual = NA:Double'])
            for ann in annotations_desc.keys():
                pvds = pvds.set_va_attributes('va.info.{}'.format(ann), annotations_desc[ann])

            pvds.export_vcf("{}.chr{}.{}.vcf.bgz".format(outprefix, interval, pname))

try_slack(['@laurent'], get_ab_stats, "gs://gnomad-lfran/tmp/gnomad.abstats")


def vep_hardcalls():
    # vds = hc.read(hardcalls_exomes_vds_path(split=True))
    # vep = hc.read("gs://gnomad-exomes-raw-hail01/gnomad.exomes.vep.split.vds")
    # vds=vds.annotate_variants_vds(vep, expr='va.vep = vds.vep')
    # vds.write(hardcalls_exomes_vds_path(split=True) + ".vep.vds")

    vds = hc.read(hardcalls_genomes_vds_path(split=True))
    vep = hc.read("gs://gnomad-genomes-raw/gnomad.genomes.vep.split.vds")
    vds = vds.annotate_variants_vds(vep, expr='va.vep = vds.vep')
    vds.write(hardcalls_genomes_vds_path(split=True) + ".vep.vds")

    vds = hc.read(hardcalls_exomes_vds_path())
    vds = vds.vep(vep_config)
    vds.write(hardcalls_exomes_vds_path() + ".vep.vds")

    vds = hc.read(hardcalls_genomes_vds_path())
    vds = vds.vep(vep_config)
    vds.write(hardcalls_genomes_vds_path() + ".vep.vds")

    vds = hc.read(hardcalls_genomes_vds_path(adj=True))
    vep = hc.read(hardcalls_genomes_vds_path() + ".vep.vds")
    vds = vds.annotate_variants_vds(vep, expr='va.vep=vds.vep')
    vds.write(hardcalls_genomes_vds_path(adj=True) + ".vep.vds")

#try_slack(['@laurent','@konradjk'], vep_hardcalls)

def write_indels():
    vds = hc.read(final_genome_split_vds_path)
    vds = vds.filter_variants_expr('v.altAllele.isIndel')
    vds.export_variants('gs://gnomad-lfran/indels/gnomad_indels.tsv.gz',
                        ",".join([
                            'chrom = v.contig',
                            'pos = v.start',
                            'ref = v.ref',
                            'alt = v.alt',
                            'wasFiltered = va.filters.isEmpty()',
                            'filters = va.filters.mkString(",")',
                            'wasSplit = va.wasSplit',
                            'AC = va.info.AC',
                            'AN = va.info.AN',
                            'most_severe_consequence = va.vep.most_severe_consequence'
                        ]))
#write_indels()

def fix_trios():

    trios = hc.read("gs://gnomad/compound_hets/exome_trios_clean.vds")
    trios = trios.annotate_variants_expr(['va = drop(va, pass)'])
    trios = trios.annotate_variants_expr(['va.filters = va.release.filters',
                                          'va.release = select(va.release, info)'])
    trios.write("gs://gnomad/compound_hets/exome_trios_clean.fixed.vds", overwrite=True)
#fix_trios()

def get_trios_pbt_kt():
    trios = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds_path))
    trios = trios.annotate_variants_vds(hc.read(final_exome_split_vds_path), root='va.release')
    ped = Pedigree.read(exomes_fam_path)
    trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')

    ped = ped.filter_to(trios.sample_ids)
    logger.info("Found {0} trios in VDS.".format(len(ped.complete_trios()) / 3))

    trios = trios.filter_variants_expr(
        'va.release.filters.isEmpty() && va.calldata.all_samples_raw.AF <= {0} && gs.filter(g => g.isCalledNonRef).count() > 0'.format(
            args.max_af))

def create_public_pca_vds():
    vds = hc.read(gnomad_pca_vds_path)
    vds = vds.drop_samples()
    vds = vds.annotate_variants_expr([
        'va = select(va,pca_loadings)'
    ])
    vds = vds.repartition(500, shuffle=True)
    pprint(vds.variant_schema)

    vds.write('gs://gnomad-public/release/2.0.1/pca/gnomad_pca_loadings.vds', overwrite=True)
#create_public_pca_vds()

def export_avg_dp():
    vds = hc.read(full_exome_01_vds_path)
    vds = vds.annotate_samples_expr("sa.mean_dp = gs.map(g => g.dp).stats().mean")
    vds.export_samples('gs://gnomad-lfran/exomes_avg_dp.txt','s = s, mean_dp = sa.mean_dp')
#export_avg_dp()

def extract_wgs_syndip_vcf():
    vds = hc.read(full_genome_vds_path)
    vds = vds.filter_samples_list(['CHMI_CHMI3_WGS1'])

    print("")

    vds = vds.annotate_variants_expr('va.keep = gs.filter(g => g.isCalledNonRef).count() > 0')
    vds = vds.filter_variants_expr('va.keep')
    vds.export_vcf('gs://gnomad-genomes/subsets/syndip/syndip.vcf.bgz')

#extract_wgs_syndip_vcf()

def load_annotated_gnomad(exomes, release_only = False, split = False, vep = False):

    if split and vep:
        exit("Currently cannot VEP nonsplit")

    if exomes:
        data = add_exomes_sa(hc.read(full_exome_vds_path))
        if release_only:
            data = data.filter_samples_expr('sa.meta.drop_status == "keep"')




def missense_in_finns(min_AC, exomes = True):

    if exomes:
        vds = hc.read(final_exome_split_vds_path)
    else:
        vds = hc.read(final_genome_split_vds_path)

    vds = process_consequences(vds)
    vds = vds.filter_variants_expr('va.vep.worst_csq == "missense_variant"')
    print(
        vds.query_variants([
            'variants.filter(v => va.info.AC_FIN >= {} && va.info.AC_FIN == va.info.AC).count()'.format(min_AC),
            #'variants.filter(v => va.info.AC_FIN >0 && va.info.AC >=2).count()',
            'variants.filter(v => va.info.AC_FIN >=3).count()'.format(min_AC)
        ])
     )
#missense_in_finns(3, True)
#missense_in_finns(3, False)


def get_brca12_comt_table(exomes = True):
    rs165631 = Variant(contig="22", start=19951816, ref="C", alts = "T")

    if exomes:
        data = add_exomes_sa(hc.read(full_exome_vds_path))
        data = data.filter_samples_expr('sa.meta.drop_status == "keep"')
    else:
        data = add_genomes_sa(hc.read(full_genome_vds_path))
        data = data.filter_samples_expr('sa.meta.keep')
        data = data.annotate_samples_expr(['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop'])

    data = data.filter_genotypes(ADJ_CRITERIA)

    brca12_vds = data.filter_intervals([Interval.parse("17:41196312-41277500 "), Interval.parse("13:32889611-32973805")])
    brca12_vds = brca12_vds.split_multi()
    brca12_vds = brca12_vds.annotate_variants_vds(hc.read(full_exomes_vep_split_vds_path), expr="va.vep = vds.vep")
    brca12_vds = process_consequences(brca12_vds)
    brca12_vds = brca12_vds.filter_variants_expr('va.vep.worst_csq_suffix[-2:] == "HC"')
    brca12_kt = brca12_vds.genotypes_table()
    brca12_kt = brca12_kt.filter("g.isCalledNonRef")
    brca12_kt = brca12_kt.aggregate_by_key("s = s", ['nBRCA12_LoF = v.count()',
                                                     'BRCA12_vars = v.map(v => str(v)).collect().mkString(",")'])

    rs165631_kt = data.filter_variants_list([rs165631]).genotypes_table()
    rs165631_kt = rs165631_kt.filter("g.isCalled")
    rs165631_kt = rs165631_kt.annotate(["pop = sa.meta.population", "rs165631 = g.isCalledNonRef"])
    rs165631_kt = rs165631_kt.select(["s","pop","rs165631"]).key_by("s")

    kt = rs165631_kt.join(brca12_kt, how="left")
    kt.export("gs://gnomad-lfran/brca12_mguo/gnomad_{}_brca12_lof_carriers.tsv".format("exomes" if exomes else "genomes"))

#get_brca12_comt_table(False)


def export_pca():
    vds = hc.read(gnomad_pca_vds_path)
    vds.export_samples("gs://gnomad-genomes/sampleqc/gnomad.pca.txt.gz",
                       'data = s.split("_")[0],underscored_sample = s.split("_")[1:].mkString("_"), pc1 = sa.pca.PC1,pc2 = sa.pca.PC2,pc3 = sa.pca.PC3,pc4 = sa.pca.PC4,'
                       'pc5 = sa.pca.PC5,pc6 = sa.pca.PC6,pc7 = sa.pca.PC7,pc8 = sa.pca.PC8,pc9 = sa.pca.PC9,'
                       'pc10 = sa.pca.PC10,pc11 = sa.pca.PC11,pc12 = sa.pca.PC12,pc13 = sa.pca.PC13,pc14 = sa.pca.PC14,'
                       'pc15 = sa.pca.PC15,pc16 = sa.pca.PC16,pc17 = sa.pca.PC17,pc18 = sa.pca.PC18,pc19 = sa.pca.PC19,'
                       'pc20 = sa.pca.PC20'
                       )
#export_pca()

def get_samples_at_variant_pairs2(variant_pair_list_path, out_path, parallel=False, sample_list=None, export_homref=False):

    def get_intervals_list(kt, v_index):
        ll = kt.annotate('l = v{}.locus()'.format(v_index)).select(["l"]).collect()
        return([Interval(x.l, Locus(x.l.contig, x.l.position + 1)) for x in ll])

    vp_kt = hc.import_table(variant_pair_list_path, impute=True)
    intervals = get_intervals_list(vp_kt,"1") + get_intervals_list(vp_kt,"2")

    vds = hc.read(full_exome_01_vds_path)
    vds = vds.filter_intervals(intervals)

    if sample_list:
        vds = vds.filter_samples_list(sample_list)

    vds = vds.split_multi()
    geno_kt = vds.genotypes_table().drop(['va','sa'])

    vp_kt = vp_kt.key_by('v1').join(geno_kt.key_by(['v']), how = "inner")
    vp_kt = vp_kt.rename({'g' : 'g1'})
    vp_kt = vp_kt.key_by(['v2', 's']).join(geno_kt.key_by(['v','s']), how="inner")
    vp_kt = vp_kt.rename({'g': 'g2'})

    if not export_homref:
        vp_kt = vp_kt.filter("g1.isCalledNonRef || g2.isCalledNonRef")

    vp_kt = vp_kt.annotate([
        'v1_gt = g1.gt',
        'v1_gq = g1.gq',
        'v1_dp = g1.dp',
        'v1_het = g1.isHet',
        'v1_ad0 = g1.ad[0]',
        'v1_ad1 = g1.ad[1]',
        'v1_pl0 = g1.pl[0]',
        'v1_pl1 = g1.pl[1]',
        'v1_pl2 = g1.pl[2]',
        'v2_gt = g2.gt',
        'v2_gq = g2.gq',
        'v2_dp = g2.dp',
        'v2_het = g2.isHet',
        'v2_ad0 = g2.ad[0]',
        'v2_ad1 = g2.ad[1]',
        'v2_pl0 = g2.pl[0]',
        'v2_pl1 = g2.pl[1]',
        'v2_pl2 = g2.pl[2]',
        'nVars = g1.isCalledNonRef.toInt() + g2.isCalledNonRef.toInt()'
    ])
    vp_kt = vp_kt.drop(['g1','g2'])
    vp_kt.export(out_path, parallel=parallel)

#hadoop_copy('gs://gnomad-exomes-hail01/variantqc/exac2.qctrios.fam', 'file:///exac2.qctrios.fam')
# get_samples_at_variant_pairs2("gs://gnomad/compound_hets/gnomad_no_pop.discordance_phase_vp.tsv",
#                               "gs://gnomad/compound_hets/gnomad_no_pop.discordance_phase_vp.exomes_trios_samples.tsv",
#                               sample_list=[trio.proband for trio in Pedigree.read("gs://gnomad-exomes-hail01/variantqc/exac2.qctrios.fam").trios],
#                               export_homref=True)


def get_non_ref_samples_at_variant_pairs(variants, out):
    vds = hc.read(full_exome_vds_path)
    vds = vds.filter_variants_list(variants)
    vds = vds.annotate_samples_expr("sa.nVars = gs.filter(g => g.isCalledNonRef).count()")
    vds = vds.filter_samples_expr("sa.nVars > 0")
    kt = vds.genotypes_table()
    kt = kt.flatten()
    kt = kt.annotate(
        ['gt = g.gt',
         'gq = g.gq',
         'dp = g.dp',
         'het = g.isHet',
         'ad0 = g.ad[0]',
         'ad1 = g.ad[1]',
         'pl0 = g.pl[0]',
         'pl1 = g.pl[1]',
         'pl2 = g.pl[2]']
    )
    kt.export(out)

#variants = [Variant.parse("2:179410666:G:A"), Variant.parse("2:179453458:G:A")]
#variants = [Variant.parse("16:27475815:C:T"), Variant.parse("16:27506594:G:A")]
#variants = [Variant.parse("7:100678935:C:T"), Variant.parse("7:100679429:T:A")]
#get_non_ref_samples_at_variant_pairs(variants,
#                                     'gs://gnomad-lfran/tmp/{}.tsv'.format("_".join([str(v) for v in variants])))

def get_gtHaplotypeCounts(variants, out_path):
    import compound_hets_utils

    def get_kt(vds, variants, exomes):

        variant_pairs = []
        for i, v1 in enumerate(variants[:-1]):
            for v2 in variants[i+1:]:
                variant_pairs.append({"v1": v1, "v2": v2})

        variants_kt = KeyTable.from_py(hc,
                                       variant_pairs,
                                       TStruct(["v1","v2"], [TVariant(), TVariant()]),
                                       key_names=["v1","v2"],
                                       num_partitions=1)

        out_prefix = "ge_" if exomes else "gg_"

        if exomes:
            vds = add_exomes_sa(vds)
            vds = vds.annotate_samples_expr(['sa.pop = sa.meta.population'])
            vds = vds.filter_samples_expr('sa.meta.drop_status == "keep"')
        else:
            vds = add_genomes_sa(vds)
            vds = vds.annotate_samples_expr(['sa.pop = sa.meta.final_pop'])
            vds = vds.filter_samples_expr('sa.meta.keep')

        vds = vds.filter_intervals([Interval.parse("{}:{}-{}".format(v.contig, v.start, v.start+1)) for v in variants])
        vds = vds.split_multi()
        vds = vds.filter_variants_list(variants)

        vds = filter_to_adj(vds)
        vds = vds.annotate_variants_expr('va.contig = v.contig')
        kt = vds.phase_em(['va.contig'], num_partitions=1, sa_keys=['sa.pop'], variant_pairs=variants_kt)
        kt = kt.aggregate_by_key(key_expr=["v1=v1", "v2=v2", "pop=`sa.pop`"],
                                 agg_expr=[out_prefix + expr for expr in
                                           ["genotype_counts = genotype_counts.take(1)[0]",
                                           "haplotype_counts = haplotype_counts.take(1)[0]",
                                           "prob_same_haplotype = prob_same_haplotype.take(1)[0]"]])
        kt, _ = compound_hets_utils.flatten_haplotype_counts(kt,
                                                             gc_col=out_prefix + "genotype_counts",
                                                             hc_col=out_prefix + "haplotype_counts",
                                                             out_prefix=out_prefix)
        return kt


    exomes = get_kt(hc.read(full_exome_vds_path), variants, True)
    genomes = get_kt(hc.read(full_genome_vds_path), variants, False)

    both = exomes.join(genomes, how="outer")
    both = both.annotate(["{0} = orElse(ge_{0},0) + orElse(gg_{0},0)".format(field[3:]) for field in
                   both.columns if field.startswith("ge_") and field[3:] in compound_hets_utils.gc_cols ])
    both = both.annotate(["{0} = orElse(ge_{0},0.0) + orElse(gg_{0},0.0)".format(field[3:]) for field in
                          both.columns if field.startswith("ge_") and field[
                                                                      3:] in compound_hets_utils.hc_cols])

    both.export(out_path)

#IFIH1
    # rs35744605 => 2:163134090:C:A
    # rs35337543 => 2:163136505:C:G
    # rs35732034 => 2:163124596:C:T
    # rs35667974 => 2:163124637:T:C
IFIH1_variants = [
    Variant.parse("2:163124596:C:T"),
    Variant.parse("2:163124637:T:C"),
    Variant.parse("2:163134090:C:A"),
    Variant.parse("2:163136505:C:G")
]

#get_gtHaplotypeCounts(ifih1_variants, "gs://gnomad/compound_hets/variant_lookups/gnomad_IFIH1_LoF.counts.txt")

LEPR_variants = [
    Variant.parse("1:66081791:C:T"),
    Variant.parse("1:66102679:T:C")
]
#get_gtHaplotypeCounts(LEPR_variants, "gs://gnomad/compound_hets/variant_lookups/gnomad_LEPR.counts.txt")

def add_vqslod():
    rf = (hc.read("gs://gnomad-exomes/variantqc/170620_new/gnomad_exomes.annotated_for_rf.novqslod.vds", drop_samples=True)
          .annotate_variants_vds(hc.read(vqsr_vds_path).split_multi(),
                                 expr='va.info.VQSLOD = vds.info.VQSLOD,'
                                      'va.info.POSITIVE_TRAIN_SITE = vds.info.POSITIVE_TRAIN_SITE,'
                                      'va.info.NEGATIVE_TRAIN_SITE = vds.info.NEGATIVE_TRAIN_SITE'))
    rf.write("gs://gnomad-exomes/variantqc/170620_new/gnomad_exomes.annotated_for_rf.vds",overwrite=True)
#add_vqslod()

def reprocess_sites():
    exomes= hc.read('gs://gnomad-exomes/sites/vds/gnomad.exomes.r2.0.1.sites.split.vds')
    exomes = exomes.annotate_variants_expr([
        'va.filters = va.filters.filter(x => x != "RF" && x != "AC0").union(va.info.AS_FilterStatus)',
        'va.info.AS_RF_POSITIVE_TRAIN = va.info.AS_RF_POSITIVE_TRAIN.toSet.contains(va.aIndex)',
        'va.info.AS_RF_NEGATIVE_TRAIN = va.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(va.aIndex)'
    ])
    exomes.write(final_exome_split_vds_path, overwrite=True)

    genomes= hc.read('gs://gnomad-genomes/sites/vds/gnomad.genomes.r2.0.1.sites.split.vds')
    genomes = genomes.annotate_variants_expr([
        'va.filters = va.filters.filter(x => x != "RF" && x != "AC0").union(va.info.AS_FilterStatus)',
        'va.info.AS_RF_POSITIVE_TRAIN = va.info.AS_RF_POSITIVE_TRAIN.toSet.contains(va.aIndex)',
        'va.info.AS_RF_NEGATIVE_TRAIN = va.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(va.aIndex)'
    ])
    genomes.write(final_genome_split_vds_path, overwrite=True)

#reprocess_sites()


def tj_variant():
    path = 'gs://gnomad-exomes/sites/vds/gnomad.exomes.r2.0.1.sites.split.vds'
    path = final_exome_split_vds_path
    vds = hc.read(path)
    query = vds.filter_intervals(Interval.parse('16:30977316-30977317')).query_variants('variants.map(v => {v:v, Filter: va.filters, AS_Filter: va.info.AS_FilterStatus}).collect()')
    pprint(query)

#tj_variant()


def get_psych_sample_meta():
    vds = hc.read("gs://gnomad-sczmeta/genomes/sczmeta.vds", drop_variants=True)
    vds = add_genomes_sa(vds)
    kt = vds.samples_table()
    pprint(kt.schema)
    kt.write("gs://gnomad-sczmeta/genomes/sczmeta.samples_meta.kt")
#get_psych_sample_meta()

def merge_ann_rf(genomes_rf_ann, exomes_rf_ann):
    genomes_rf_ann = genomes_rf_ann.annotate_variants_expr('va.stats = drop(va.stats, all_samples_raw, release_samples_raw)')
    genomes_kt = genomes_rf_ann.variants_table().rename({'va': 'genomes'})
    exomes_rf_ann = exomes_rf_ann.annotate_variants_expr(
        'va.stats = drop(va.stats, all_samples_raw, release_samples_raw)')
    exomes_kt = exomes_rf_ann.variants_table().rename({'va': 'exomes'})
    merged_kt = genomes_kt.join(exomes_kt, how="outer")
    VariantDataset.from_table(merged_kt).write("gs://gnomad-genomes/variantqc/merged_exomes_genomes_for_rf.vds", overwrite=True)
#merge_ann_rf(hc.read("gs://gnomad-genomes/variantqc/rf_new_stats.annotated_for_rf.vds"),
#             hc.read("gs://gnomad-exomes/variantqc/170620_new/gnomad_exomes.annotated_for_rf.vds"))

def get_non_ref_gts():
    genomes = hc.read(full_genome_vds_path).filter_samples_list(
        read_list_data(genomes_qc_pass_samples_list_path))
    genomes = genomes.filter_intervals(Interval.parse("1:871215-871216"))
    pprint(genomes.query_genotypes("gs.filter(g => g.isCalledNonRef).map(g => {s:s, g:g}).collect()"))
#get_non_ref_gts()

def check_variant():
    vds = hc.read(full_exome_vds_path)
    chrom = 5
    pos1 = 118513128
    pos2 = 118573017
    #vds = vds.filter_samples_expr('sa.meta.drop_status == "keep"')
    vds = vds.filter_samples_list(['MH0132077', 'MH0132083', 'MH0132078'])
    vds = vds.filter_intervals([Interval.parse(x) for x in ['{}:{}-{}'.format(chrom, pos1, pos1+1),
                                                            '{}:{}-{}'.format(chrom, pos2, pos2+1)]])
    #vds = vds.filter_samples_list(['MH0126636','MH0126637','MH0126638'])
    #vds = vds.filter_intervals([Interval.parse(x) for x in ['9:139945737-139945738','2:139946791-139946792']])
    pprint(vds.query_genotypes(['gs.map(g => {v:v, s:s, g:g}).collect()']))

#check_variant()

def get_bravo_fingerprint():

    print("Computing gnomAD fingerprint")

    #sites = hc.import_vcf("gs://gnomad-resources/gnomad-bravo/ExAC.r1.500.fingerprint.sites.vcf.gz", force=True)
    #sites.write("gs://gnomad-resources/gnomad-bravo/ExAC.r1.500.fingerprint.sites.vds")
    sites = hc.read("gs://gnomad-resources/gnomad-bravo/ExAC.r1.500.fingerprint.sites.vds").variants_table().select(["v"])
    exomes = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds_path)).filter_variants_table(sites)
    print("Found {} samples in exomes before filtering.".format(len(exomes.sample_ids)))
    exomes = exomes.filter_samples_expr('sa.meta.drop_status == "keep"')
    exome_samples = exomes.sample_ids
    print("Found {} samples in exomes after filtering.".format(len(exome_samples)))
    exomes = exomes.rename_samples({old: "sample_{}".format(new) for new,old in enumerate(exome_samples)})
    exomes = exomes.annotate_variants_expr('va = {}')
    exomes = exomes.annotate_samples_expr('sa = {}')
    exomes.persist()
    print("Exome variants overlapping fingerprint sites: {}".format(exomes.count_variants()))
    genomes = add_genomes_sa(hc.read(full_genome_hardcalls_split_vds_path).filter_variants_table(sites))
    print("Found {} samples in genomes before filtering.".format(len(genomes.sample_ids)))
    genomes = genomes.filter_samples_expr('sa.meta.keep')
    print("Found {} samples in genomes after filtering.".format(len(genomes.sample_ids)))
    genomes = genomes.rename_samples(
        {old: "sample_{}".format(new + len(exome_samples)) for new, old in enumerate(genomes.sample_ids)})
    genomes = genomes.annotate_variants_expr('va = {}')
    genomes = genomes.annotate_samples_expr('sa = {}')
    genomes.persist()
    print("Genome variants overlapping fingerprint sites: {}".format(genomes.count_variants()))
    both = exomes.join(genomes)
    both.persist()
    print("Combined variants overlapping fingerprint sites: {}".format(both.count_variants()))
    both.export_vcf("gs://gnomad-resources/gnomad-bravo/gnomad_samples.fingerprint.vcf.bgz")

#get_bravo_fingerprint()




#print(hc.grep('^22\t165561[1246][92]\t','gs://gnomad-genomes/subsets/afib/afib.22.vcf.bgz'))

def create_release_split_vds(vds, out):
    annotations, a_annotations, g_annotations, dot_annotations = get_numbered_annotations(vds, "va.info")

    vds = vds.split_multi()
    vds = vds.annotate_variants_expr(index_into_arrays(a_based_annotations=["va.info." + a.name for a in a_annotations], vep_root='va.vep') )
    vds = vds.annotate_variants_expr(['va.filters = va.filters.filter(x => !["AC0","RF"].toSet.difference(va.info.AS_FilterStatus).contains(x))',
                                     'va.info.AS_RF_NEGATIVE_TRAIN = isDefined(va.info.AS_RF_NEGATIVE_TRAIN) && va.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(va.aIndex)',
                                      'va.info.AS_RF_POSITIVE_TRAIN = isDefined(va.info.AS_RF_POSITIVE_TRAIN) && va.info.AS_RF_POSITIVE_TRAIN.toSet.contains(va.aIndex)',
                                      'va.info = drop(va.info, {0} )'.format(",".join([a.name for a in g_annotations]))])

    vds.write(out)

#create_release_split_vds(hc.read(final_exome_vds), final_exome_split_vds)
#create_release_split_vds(hc.read(final_genome_vds), final_genome_split_vds)




def rename_samples_raw():
    vds = hc.read('gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.new_stats.vds')
    vds = vds.annotate_variants_expr([
        'va.calldata.all_samples_raw = va.calldata.raw'

    ])
    vds = vds.annotate_variants_expr([
        'va.calldata = drop(va.calldata, raw)'

    ])
    vds.write(full_genome_hardcalls_vds_path, overwrite = True)

#rename_samples_raw()

def recompute_qc_samples_stats():

    sample_group_filters = {"qc_samples_raw": 'sa.meta.qc_sample || (sa.in_exomes && sa.qc_pass)'}

    allele_annotations = []
    for group, filter_expr in sample_group_filters.iteritems():
        allele_annotations.extend(
            get_allele_stats_expr("va.stats.%s" % group, medians=True, samples_filter_expr=filter_expr))

    variant_annotations = []
    for group, filter_expr in sample_group_filters.iteritems():
        if filter_expr:
            filter_expr = '.filter(g => %s)' % filter_expr

        variant_annotations.extend(["va.calldata.%s = gs%s.callStats(g => v)" % (group, filter_expr),
                                    "va.stats.%s.qd = (va.stats.%s.qual / va.stats.%s.nrdp).map(x => [35,x].min)" % (
                                        group, group, group)])

    vds = add_genomes_sa(hc.read(full_genome_vds_path))
    vds = vds.annotate_variants_vds(hc.read(kgp_high_conf_snvs_vds_path), 'va.kgp_high_conf = isDefined(vds)')
    (
        vds.annotate_alleles_expr(allele_annotations)
            .annotate_variants_expr(variant_annotations)
            .drop_samples()
            .write("gs://gnomad-genomes/tmp/new_qc_samples_annotations.vds", overwrite=True)
    )

    annotations = hc.read("gs://gnomad-genomes/tmp/new_qc_samples_annotations.vds")
    vds = hc.read(full_genome_hardcalls_vds_path)
    vds.annotate_variants_vds(annotations, 'va.calldata.qc_samples_raw = vds.calldata.qc_samples_raw,'
                                           'va.stats.qc_samples_raw = vds.stats.qc_samples_raw')
    vds.write("gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.new_qc_samples_stats.vds")

    # vds = hc.read(full_genome_hardcalls_split_vds)
    # vds.annotate_variants_vds(annotations.split_multi(), 'va.calldata.qc_samples_raw = vds.calldata.qc_samples_raw,'
    #                                                      'va.stats.qc_samples_raw = vds.stats.qc_samples_raw')
    # vds.write("gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.new_qc_samples_stats.split.vds")
#recompute_qc_samples_stats()


def export_clinvar_gts(out = "gs://gnomad-lfran/tmp/clinvar.gg_filtered.genotypes.txt.bgz"):
    gg = hc.read(full_genome_vds_path)
    variants = [hail.Variant.parse(v) for v in read_list_data("gs://gnomad-lfran/tmp/clinvar.gg_filtered.txt")[1:]]
    gg = gg.filter_variants_list(variants).filter_multi()
    gg = gg.annotate_variants_expr('va.nNonRef = gs.filter(g => g.isCalledNonRef).count(), '
                              'va.nMissing = gs.filter(g => !g.isCalled).count()')
    gg.export_genotypes(out, ",".join([
        'chrom = v.contig',
        'pos = v.start',
        'ref = v.ref',
        'alt = v.alt',
        'sample = s',
        'gt = g.gt',
        'gq = g.gq',
        'dp = g.dp',
        'pl0 = g.pl[0]',
        'pl1 = g.pl[1]',
        'pl2 = g.pl[2]',
        'pab = g.pAB()',
        'ad0 = g.ad[0]',
        'ad1 = g.ad[1]',
        'nNonRef = va.nNonRef',
        'nMissing = va.nMissing'
    ]))
#export_clinvar_gts()

def correct_qual_for_rf():
    rf = hc.read("gs://gnomad-genomes/variantqc/rf_as_qd_only.annotated_for_rf.vds")
    qd = hc.read(full_genome_hardcalls_split_vds_path)

    sample_group_filters = {"all_samples_raw": '',
                            "qc_samples_raw": 'sa.meta.qc_sample',
                            "release_samples_raw": 'sa.meta.keep'
                            }
    variant_expr = []
    for group in sample_group_filters.keys():
        variant_expr.extend([
            'va.stats.%s.qd = vds.stats.%s.qd[vds.aIndex - 1]' % (group, group),
            'va.stats.%s.qual = vds.stats.%s.qual[vds.aIndex - 1]' % (group, group)
        ])


    rf = rf.annotate_variants_vds(qd, expr = ",".join(variant_expr))
    rf.write("gs://gnomad-genomes/variantqc/rf_new_stats.annotated_for_rf.vds")

#correct_qual_for_rf()

def compute_qd(vds, out):
    sample_group_filters = {"all_samples_raw": '',
                            "qc_samples_raw": 'sa.meta.qc_sample',
                            "release_samples_raw": 'sa.meta.keep'
                            }

    variant_expr = []
    for group in sample_group_filters.keys():
        variant_expr.extend([
            'va.stats.%s.qd = (va.stats.%s.qual / va.stats.%s.nrdp).map(x => [x,35].min)' % (group, group, group)
        ])

    vds.annotate_variants_expr(variant_expr).write(out, overwrite = True)

#compute_qd(hc.read("gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.noqd.vds"), full_genome_hardcalls_vds)
#compute_qd(hc.read("gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.split.noqd.vds"), full_genome_hardcalls_split_vds)

def add_new_features_to_hardcalls(exomes = False):

    logger.info("Adding new features to Hardcalls")

    def get_new_allele_stats_expr(root="va.stats", samples_filter_expr=''):

        if samples_filter_expr:
            samples_filter_expr = "&& " + samples_filter_expr

        stats = [
                 '%s.pab = gs.filter(g => g.isHet %s).map(g => g.pAB()).stats()',
                 '%s.nrdp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).sum()',
                 '%s.qual = -10*gs.filter(g => g.isCalledNonRef %s).map(g => if(g.pl[0] > 3000) -300 else log10(g.gp[0])).sum()',
                 '%s.combined_pAB = let hetSamples = gs.filter(g => g.isHet %s).map(g => log(g.pAB())).collect() in orMissing(!hetSamples.isEmpty, -10*log10(pchisqtail(-2*hetSamples.sum(),2*hetSamples.length)))',
                 '%s.pab_median = gs.filter(g => g.isHet %s).map(g => g.pAB()).collect().median']

        stats_expr = [x % (root, samples_filter_expr) for x in stats]

        return stats_expr

    def create_new_annotations(vds, sample_group_filters):
        allele_annotations = []
        for group, filter_expr in sample_group_filters.iteritems():
            allele_annotations.extend(
                get_new_allele_stats_expr("va.stats.%s" % group, samples_filter_expr=filter_expr))

        variant_annotations = []
        for group, filter_expr in sample_group_filters.iteritems():
            variant_annotations.extend(["va.stats.%s.qd = (va.stats.%s.qual / va.stats.%s.nrdp).map(x => [35,x].min)" % (group, group, group)])

        hapmap = hc.read(hapmap_vds_path)
        mills = hc.read(mills_vds_path)
        omni = hc.read(omni_vds_path)

        (
            vds
                .annotate_variants_vds(hapmap, expr='va.hapmap = isDefined(vds)')
                .annotate_variants_vds(omni, expr='va.omni = isDefined(vds)')
                .annotate_variants_vds(mills, expr='va.mills = isDefined(vds)')
                .annotate_alleles_expr(allele_annotations)
                .annotate_variants_expr(variant_annotations)
                .drop_samples()
                .write("gs://gnomad-lfran/genomad_genomes.hardcalls_new_ann.vds", overwrite=True)
        )

    def annotate_with_new_annotations(hardcalls, new_annotations,  sample_group_filters, out):

        hardcalls.annotate_variants_expr(['va.calldata.all_samples_raw = va.calldata.raw',
                                          'va.calldata = drop(va.calldata, raw)'])

        vds_annotations = ['va.hapmap = vds.hapmap', 'va.omni = vds.omni', 'va.mills = vds.mills']
        for group in sample_group_filters.keys():
            for metric in ['pab', 'nrdp', 'qual', 'combined_pAB', 'pab_median']:
                vds_annotations.append('va.stats.%s.%s = vds.stats.%s.%s' % (group, metric, group, metric))

        hardcalls = hardcalls.annotate_variants_vds(new_annotations, ",".join(vds_annotations))
        hardcalls.write(out, overwrite=True)

    sample_group_filters = {"all_samples_raw": '',
                            "qc_samples_raw": 'sa.meta.qc_sample',
                            "release_samples_raw": 'sa.meta.keep'
                            }

    if exomes:
        vds = (
            hc.read(full_genome_vds_path)
                .annotate_samples_table(hail.KeyTable.import_fam(genomes_fam_path), root='sa.fam')
                .annotate_samples_table(hc.import_table(genomes_meta_tsv_path, impute=True).key_by('Sample'), root='sa.meta')
        )
    else:
        vds = (
            hc.read(full_genome_vds_path)
                .annotate_samples_table(hail.KeyTable.import_fam(genomes_fam_path), root='sa.fam')
                .annotate_samples_table(hc.import_table(genomes_meta_tsv_path, impute=True).key_by('Sample'), root='sa.meta')
        )

    create_new_annotations(vds, sample_group_filters)
    new_annotations = hc.read("gs://gnomad-lfran/genomad_genomes.hardcalls_new_ann.vds")
    annotate_with_new_annotations(hc.read("gs://gnomad-lfran/tmp/gnomad.raw_hardcalls.vds"),
                                  new_annotations,
                                  sample_group_filters,
                                  full_genome_hardcalls_vds_path)

    annotate_with_new_annotations(hc.read("gs://gnomad-lfran/tmp/gnomad.raw_hardcalls.split.vds"),
                                  new_annotations.split_multi(),
                                  sample_group_filters,
                                  full_genome_hardcalls_split_vds_path)


#add_new_features_to_hardcalls()

def find_dups():
    logger.info("Finding duplicate samples")

    genomes = hc.read(full_genome_hardcalls_vds_path)
    exomes = hc.read(full_exome_vds_path)

    samples = list(set(genomes.sample_ids).intersection(set(exomes.sample_ids)))
    logger.debug("Found %d samples with same ID." % len(samples))

    genomes = genomes.filter_samples_list(samples).rename_samples({x: "gg_"+x for x in samples})
    exomes = exomes.filter_samples_list(samples).rename_samples({x: "ge_" + x for x in samples})

    p5k = hail.KeyTable.import_interval_list(purcell5k_intervals_path)
    genomes = genomes.filter_variants_table(p5k).filter_multi()
    exomes = exomes.filter_variants_table(p5k).filter_multi()

    genomes = genomes.annotate_samples_expr('sa = {}')
    exomes = exomes.annotate_samples_expr('sa = {}')

    combined = exomes.join(genomes)

    ibd = combined.ibd(maf = 'va.info.AF[0]', min = 0.7)
    ibd.write('gs://gnomad-lfran/tmp/gnomad-ibd.kt')

    ibd.export('gs://gnomad-lfran/tmp/gnomad-ibd.txt.gz')



#
# vds = hc.read(full_genome_vds)
# vds = vds.filter_samples_expr('s == "NWD108733"')
# vds = vds.filter_variants_intervals(Interval.parse("22:16556118-16556163"))
# print(vds.query_genotypes('gs.map(g => {v: v, g:g}).collect()'))


#hc.read(full_genome_vds, sites_only=True).split_multi().vep(config =vep_config).write(full_genomes_vep_split_vds, overwrite=True)
#hc.read(full_exome_vds, sites_only=True).split_multi().vep(config =vep_config).write(full_exomes_vep_split_vds, overwrite=True)

def migraine_subset():
    logger.info("Generating migraine subset for John")
    samples = read_list_data('gs://gnomad-genomes/subsets/migraine/migraine_samples.txt')
    vds = hc.read(full_genome_vds_path).filter_samples_list(samples)
    vds = vds.annotate_variants_expr('va.calldata = gs.callStats(g => v)')
    vds = vds.filter_alleles('va.calldata.AC[aIndex]>0', annotation='va.aIndices = aIndices', subset= True)
    vds.write('gs://gnomad-genomes/subsets/migraine/migraine.vds')

#migraine_subset()

def get_vqslod_thresholds():
    logger.info("Getting VQSLOD thresholds")
    vds = hc.read(full_genome_vds_path, sites_only=True)
    kt = vds.filter_variants_intervals(Interval.parse("1")).variants_keytable()
    kt = kt.flatten().select(['va.filters','va.info.VQSLOD','v'])
    kt = kt.annotate(['filter = `va.filters`.toArray','indel = v.altAlleles.exists(a => a.isIndel)'])
    kt = kt.select(['filter','va.info.VQSLOD','indel']).explode('filter').aggregate_by_key('filter = filter, indel = indel','min_VQSLOD = min(`va.info.VQSLOD`)')
    kt.to_dataframe().show()
    print(kt.collect())


# vds = hc.read(full_genome_vds).filter_variants_intervals((Interval.parse('22:16423491-16423493')))
# vds.write('gs://gnomad-lfran/tmp/1var.gnomad.vds')


#run_sites_sanity_checks(hc.read(final_genome_vds), genome_pops)

def check_variant_count(name, vds, vds2):
    x = vds.query_variants('variants.count()')
    y = vds2.query_variants('variants.count()')
    if x ==y:
        print "SUCCESS!"
    else:
        print "FAILURE"

    print "VDS1: %d, VDS2: %d" % (x,y)

# check_variant_count(
#     "Afib",
#     hc.read("gs://gnomad-genomes/subsets/afib/afib.vds", sites_only=True).filter_variants_intervals((Interval.parse('22'))),
#     hc.import_vcf("gs://gnomad-genomes/subsets/afib/afib.22.vcf.bgz", sites_only=True).filter_variants_intervals((Interval.parse('22')))
# )

# check_variant_count(
#     "Odonovan",
#     hc.read("gs://gnomad-genomes/subsets/odonovan_quads/odonovan_quads.vds", sites_only=True).filter_variants_intervals((Interval.parse('22'))),
#     hc.import_vcf("gs://gnomad-genomes/subsets/odonovan_quads/odonovan_quads.22.vcf.bgz", sites_only=True).filter_variants_intervals((Interval.parse('22')))
# )

# genomes = hc.read(full_genome_vds).filter_variants_intervals((Interval.parse('22')))
# subset_vds = hc.read("gs://gnomad-sczmeta/genomes/sczmeta.vds").filter_variants_intervals((Interval.parse('22')))
# run_samples_sanity_checks(subset_vds, genomes, verbose=True, n_samples=2)
#
# genomes = genomes.split_multi(keep_star_alleles = True)
# subset_vds = subset_vds.split_multi(keep_star_alleles = True)
# run_samples_sanity_checks(subset_vds, genomes, verbose=True, n_samples=2)
#
# ref_gt_expr = ['%s_ref = g.%s' % (x,x) for x in ['gt','ad','dp','gq']]

def get_dels(vds,name ,sampleid = "00313737"):
    vds = vds.filter_samples_expr('s == "00313737"')
    vds = vds.filter_variants_expr('gs.filter(g => s == "00313737" && g.isCalledNonRef && '
                          '( (g.gtk > 0 && v.altAlleles[g.gtk-1].isDeletion) || '
                          '(g.gtj > 0 && v.altAlleles[g.gtj-1].isDeletion) )).count() > 0')

    gt_expr = '%%s_%s = g.%%s' % name
    kt = vds.make_keytable('v_%s = v' % name, [gt_expr % (x, x) for x in ['gt', 'ad', 'dp', 'gq']])
    kt = kt.annotate(['chrom = v_%s.contig' % name, 'pos =  v_%s.start' % name])
    return kt.key_by(['chrom','pos']).to_dataframe()

    # .query_genotypes('gs.map(g => {v: v, g: g}).collect()')



def get_kt(vds, sampleid):
    vds = vds.filter_samples_expr('s == "%s"' % sampleid)
    # vds = vds.filter_variants_expr(
    #      'gs.filter(g => v.altAlleles.exists(a => a.isDeletion) && g.isCalledNonRef()).count() > 0')
    vds = vds.filter_variants_expr(
         'gs.filter(g => g.isCalledNonRef()).count() > 0')
    kt = vds.variants_keytable()
    kt = kt.annotate(['chrom = v.contig', 'pos =  v.start'])
    return kt.key_by(['chrom', 'pos']).select(['chrom', 'pos', 'v'])

def get_bad_sites(vds, ref_vds, sampleid):

    vds_df = get_kt(vds, sampleid).rename({'v': 'v_subset'}).to_dataframe()
    ref_df = get_kt(ref_vds, sampleid).rename({'v': 'v_ref'}).to_dataframe()

    full_df = vds_df.join(ref_df, how='outer', on =['chrom','pos'])
    print("In ref, not in subset:")
    #full_df.where("`v_subset.contig` is NULL").show()
    print("In subset, not in ref:")
    full_df.where("`v_ref.contig` is NULL").show()


    #kt = vds_kt.join(ref_kt, how="outer")
    # print(kt.count_rows())
    # print(kt.schema)
    # print kt.query('v_subset.collect()')
    #print kt.query('v_subset.filter(x => isMissing(v_ref)).collect()[:10]')
    #return kt

#get_bad_sites(subset_vds, genomes, '01C05110')


# for i in range(len(s[1])):
#     if str(s[1][i]) != str(s[2][i]):
#         print("Failed: %s != %s" %(str(s[1][i]), str(s[2][i])))
#     else:
#         print("Success: %s != %s" % (str(s[1][i]), str(s[2][i])))



# vds = vds.filter_variants_intervals(Interval.parse("22"))
# vds = vds.filter_alleles('v.altAlleles[aIndex -1].isStar && isMissing(va.info.AS_FilterStatus[aIndex - 1])',annotation='va.info.AS_FilterStatus = va.info.AS_FilterStatus[aIndices[1] - 1]', keep_star= True)
# print(vds.query_variants(['variants.map(v => va.filters.mkString(",")).counter()',
#                           'variants.map(v => va.info.AS_FilterStatus.mkString(",")).counter()']))

#[{u'': 13643L, u'SEGDUP,LCR,AC0': 22L, u'RF,AC0,InbreedingCoeff,LCR,SEGDUP': 4L, u'SEGDUP,LCR,AC0,RF': 538L, u'InbreedingCoeff,SEGDUP': 63L, u'SEGDUP': 3998L, u'LCR,RF,AC0': 1765L, u'InbreedingCoeff,RF,AC0': 7L, u'InbreedingCoeff,SEGDUP,LCR': 72L, u'SEGDUP,LCR,RF,AC0': 389L, u'InbreedingCoeff,LCR': 171L, u'AC0,RF': 2263L, None: 21365L, u'SEGDUP,AC0,RF': 461L, u'LCR,AC0,RF': 1667L, u'RF,AC0': 3504L, u'SEGDUP,RF,AC0': 524L, u'InbreedingCoeff,SEGDUP,AC0,RF': 6L, u'SEGDUP,AC0': 25L, u'InbreedingCoeff,LCR,RF,AC0': 2L, u'InbreedingCoeff,SEGDUP,RF,AC0': 9L, u'LCR,AC0': 78L, u'LCR': 44171L, u'InbreedingCoeff': 109L, u'SEGDUP,LCR': 4107L, u'InbreedingCoeff,SEGDUP,LCR,AC0': 1L, u'AC0': 71L}, {u'AC0': 14743L, None: 84292L}]

def get_metrics(vds, output, vds_is_rf = False):
    vds = vds.split_multi()
    if not vds_is_rf:
        vds = vds.annotate_variants_vds(hc.read(genomes_rf_vds_path), code ='va.mills = vds.mills, va.hapmap = vds.hapmap, va.omni = vds.omni, va.transmitted_singleton = vds.transmitted_singleton, va.label = vds.label')

    a_ann = ['va.info.AC']
    if not vds_is_rf:
        a_ann.extend(['va.info.DP_MEDIAN', 'va.info.DREF_MEDIAN', 'va.info.GQ_MEDIAN', 'va.info.AB_MEDIAN','va.info.AS_FilterStatus'])

    export_cols = ['chrom = v.contig', 'pos = v.start', 'ref = v.ref', 'alt = v.alt', 'sor = va.info.SOR',
                   'inbreedingcoeff = va.info.InbreedingCoeff',
                   'snv = v.altAllele.isSNP',
                   'ti = v.altAllele.isTransition',
                   'tv = v.altAllele.isTransversion',
                   'mills = va.mills',
                   'hapmap = va.hapmap',
                   'omni = va.omni',
                   'transmitted_singleton = va.transmitted_singleton',
                   'ac = va.info.AC',
                   'label = va.label'
                   ]

    if vds_is_rf:
        export_cols.extend([
            'dp_median = va.stats.qc_samples_raw.dp_median',
            'dref_median = va.stats.qc_samples_raw.nrq_median',
            'gq_median = va.stats.qc_samples_raw.gq_median',
            'ab_median = va.stats.qc_samples_raw.ab_median',
            'RF_score = va.RF1.probability["TP"]'])
    else:
        export_cols.extend([
            'dp_median = va.info.DP_MEDIAN',
            'dref_median = va.info.DREF_MEDIAN',
            'gq_median = va.info.GQ_MEDIAN',
            'ab_median = va.info.AB_MEDIAN',
            'RF_filtered = va.info.AS_FilterStatus.contains("RF")'])

    vds = vds.annotate_variants_expr(index_into_arrays(a_ann))
    vds.export_variants(output, ",".join(export_cols))
#get_metrics(hc.read(final_genome_vds).filter_variants_intervals(Interval.parse("22")), 'gs://gnomad-lfran/tmp/gnomad.chr22.annotations.txt.bgz')
#get_metrics(hc.read(genomes_rf_path).filter_variants_intervals(Interval.parse("22")), 'gs://gnomad-lfran/tmp/gnomad.chr22.rf-annotations.txt.bgz', vds_is_rf = True)

def get_titv_per_class():
    hc.read(genomes_rf_vds_path).variants_keytable()

def one_variant_test():
    raw = (
        hc.read('gs://gnomad/gnom.ad.vds')
            .filter_variants_intervals(Interval.parse("22:26213360-26213361"))
            .filter_samples_expr('s == "E00859946"')
            .min_rep()
    )

    print("Raw original variants and genotype: %s" % raw.query_genotypes('gs.map(g => {v: v, g: g}).collect()'))

    raw = (
        raw.filter_alleles('v.altAlleles[aIndex -1].isIndel()')
    )

    print("Raw after filter_alleles: %s" % raw.query_genotypes('gs.map(g => {v: v, g: g}).collect()'))

    subset = (
        hc.read("gs://gnomad-sczmeta/genomes/sczmeta.vds")
            .filter_variants_intervals(Interval.parse("22:26213360-26213361"))
            .filter_samples_expr('s == "E00859946"')
    )

    print("Subset original variants and genotype: %s" % subset.query_genotypes('gs.map(g => {v: v, g: g}).collect()'))

    subset = (
        subset.filter_alleles('v.altAlleles[aIndex -1].isIndel()')
    )

    print("Subset after filter_alleles: %s" % subset.query_genotypes('gs.map(g => {v: v, g: g}).collect()'))

    global_concordance, samples_vds, variants_vds = raw.split_multi().concordance(subset.split_multi())

    print(variants_vds.query_variants('variants.collect()'))
    print(raw.query_genotypes('gs.map(g => {v: v, g: g}).collect()'))
    print(subset.query_genotypes('gs.map(g => {v: v, g: g}).collect()'))

#one_variant_test()



def get_sample_vds(vds, sampleid, contig = None, keep_star = True):
    vds = vds.filter_samples_expr('s=="%s"' % sampleid)

    if contig:
        vds = vds.filter_variants_intervals(Interval.parse(contig))

    vds = vds.min_rep()

    if not keep_star:
        vds = vds.filter_variants_expr('!v.altAlleles.exists(a => a.isStar)')

    vds = vds.filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')

    # if indelsOnly:
    #     vds = vds.filter_alleles('v.altAlleles[aIndex -1].isIndel()')

    return vds


def get_per_sample_counts(vds, sampleid, contig = None, indelsOnly = False, split_multi = False):
    vds = get_sample_vds(vds,sampleid, contig)
    if split_multi:
        vds = vds.split_multi(keep_star_alleles=True)
    return(
            vds
            .sample_qc()
            .query_samples(['samples.map(s => sa.qc.nSNP).collect()[0]',
                            'samples.map(s => sa.qc.nInsertion).collect()[0]',
                            'samples.map(s => sa.qc.nDeletion).collect()[0]'])
    )

def get_concordance(out, indelsOnly = True, keep_star = True):

    raw = get_sample_vds(hc.read(full_genome_vds_path),
                         "E00859946", contig="22", keep_star= keep_star).split_multi()
    sczmeta = (
        get_sample_vds(hc.read("gs://gnomad-sczmeta/genomes/sczmeta.vds"),
                             "E00859946", contig="22",  keep_star= keep_star)
            .split_multi()
            .annotate_variants_expr('va = {}')
        )

    if indelsOnly:
        raw = raw.filter_variants_expr('v.altAllele.isIndel')
        sczmeta = sczmeta.filter_variants_expr('v.altAllele.isIndel')


    global_concordance, samples_vds, variants_vds  = raw.concordance(sczmeta)
    pprint(global_concordance)

    variants_vds.write(out, overwrite=True)

def get_discordant_sites(vds, out, only_non_missing = True):
    filter_expr = 'let l = range(5).find(i => va.concordance[i].exists(x => x > 0)) in va.concordance[l][l] == 0'
    if only_non_missing:
        filter_expr += ' && l > 0 && va.concordance.exists(x => range(5).find(i => x[i] > 0) > 0)'
    vds = vds.filter_variants_expr(filter_expr)
    print(vds.query_variants('variants.count()'))
    vds.export_variants(out,'v = v, concordance = va.concordance')

def get_sample_gt_locus(vds, sample, locus):
    vds = vds.filter_variants_intervals(Interval.parse(locus)).filter_samples_expr('s == "%s"' % sample).split_multi()
    print(vds.query_genotypes('gs.map(g => {v: v, g: g}).collect()[0:20]'))

# with_star_path = "gs://gnomad-lfran/tmp/E00859946.concordance.chr22.indels.vds"
# no_star_path = "gs://gnomad-lfran/tmp/E00859946.concordance.chr22.indels.nostar.vds"
#
# get_concordance(with_star_path)
# get_discordant_sites(hc.read(with_star_path), 'gs://gnomad-lfran/tmp/E00859946.chr22.discordant.indels.txt', False)
#get_discordant_sites(hc.read(with_star_path), 'gs://gnomad-lfran/tmp/E00859946.chr22.discordant.indels.nonmissing.txt')

#get_concordance(no_star_path, keep_star=True)
#get_discordant_sites(hc.read(no_star_path), 'gs://gnomad-lfran/tmp/E00859946.chr22.discordant.indels.filter_variants.txt')

#get_sample_gt_locus(hc.read(full_genome_vds), "E00859946", "22:26213360-26213361")
#get_sample_gt_locus(hc.read("gs://gnomad-sczmeta/genomes/sczmeta.vds"), "E00859946", "22:26213360-26213361")


# vds = vds.annotate_variants_expr(get_variant_type_expr('va.final_variantType') + ',va.nAltAlleles = v.nAltAlleles')
# vds = vds.filter_variants_expr('va.nAltAlleles > 1 && va.final_variantType == "snv"')
# vds = vds.min_rep()
# print("\nSNPs-only\n")
# print(vds.filter_alleles('v.altAlleles[aIndex -1].isSNP').query_variants(['variants.map(v => str(v)).collect()[0:10]',
#                           'variants.map(v => v.altAlleles.map(a =>str( a.isStar)).mkString(", ")).collect()[0:10]',
#                           'variants.map(v => v.altAlleles.map(a =>str( a.isSNP)).mkString(", ")).collect()[0:10]']))
# print("\nStar-only, don't keep\n")
# print(vds.filter_alleles('v.altAlleles[aIndex -1].isStar', keep_star=False).query_variants(['variants.map(v => str(v)).collect()[0:10]',
#                           'variants.map(v => v.altAlleles.map(a =>str( a.isStar)).mkString(", ")).collect()[0:10]',
#                           'variants.map(v => v.altAlleles.map(a =>str( a.isSNP)).mkString(", ")).collect()[0:10]']))
# print("\nStar-only,keep\n")
# print(vds.filter_alleles('v.altAlleles[aIndex -1].isStar', keep_star=True).query_variants(['variants.map(v => str(v)).collect()[0:10]',
#                           'variants.map(v => v.altAlleles.map(a =>str( a.isStar)).mkString(", ")).collect()[0:10]',
#                           'variants.map(v => v.altAlleles.map(a =>str( a.isSNP)).mkString(", ")).collect()[0:10]']))
# print(vds.query_variants(['variants.map(v => v.altAlleles.map(a => a.alt).mkString(", ")).counter()',
#                     'variants.map(v => va.filters.mkString(", ")).counter()']))


#
# x = hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.X.vds")
# y = hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.X.vds")
# xy = x.union([y],merge_schemas=True)
#
# sanity_check_text = run_sanity_checks(xy.filter_variants_intervals(Interval.parse("X")),pops = POPS, contig = 'X', return_string=True, skip_star=True)
# print(sanity_check_text)
# sanity_check_text = run_sanity_checks(xy.filter_variants_intervals(Interval.parse("Y")),pops = POPS, contig='Y', return_string=True, skip_star=True)
# print(sanity_check_text)

# exomes = hc.read(final_exome_vds)
# sanity_check_text = run_sanity_checks(exomes.filter_variants_intervals(IntervalTree.parse_all(['1-22'])),pops = POPS, return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'exomes merged file -- autosomes sanity check')
# sanity_check_text = run_sanity_checks(exomes.filter_variants_intervals(Interval.parse("X")),pops = POPS, contig = 'X', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'exomes merged file -- X sanity check')
# sanity_check_text = run_sanity_checks(exomes.filter_variants_intervals(Interval.parse("Y")),pops = POPS, contig='Y', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'exomes merged file -- Y sanity check')

# genomes = hc.read(final_genome_vds)
# sanity_check_text = run_sanity_checks(genomes.filter_variants_intervals(IntervalTree.parse_all(['1-22'])),pops = genome_pops, return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'genomes merged file -- autosomes sanity check')
# sanity_check_text = run_sanity_checks(genomes.filter_variants_intervals(Interval.parse("X")),pops = genome_pops, contig = 'X', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'genomes merged file -- X sanity check')
# sanity_check_text = run_sanity_checks(genomes.filter_variants_intervals(Interval.parse("Y")),pops = genome_pops, contig='Y', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'genomes merged file -- Y sanity check')

# exomes_vdses = [
# hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.autosomes.vds"),
# hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.X.vds"),
# hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.Y.vds")
# ]
# exomes_vdses = merge_schemas(exomes_vdses)
# exomes_vdses[0].union(exomes_vdses[1:]).write(final_exome_vds)
#
# genomes_vdses = [
# hc.read("gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.autosomes.vds"),
# hc.read("gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.X.vds")
# ]
# genomes_vdses = merge_schemas(genomes_vdses)
# genomes_vdses[0].union(genomes_vdses[1:]).write(final_genome_vds)


#
# print(vds.variant_schema)
# print(vds.sample_schema)
# print(vds.query_samples(['samples.count()']))
# print(vds.query_variants(['variants.map(v => va.calldata.all_samples_raw.AN).stats().max',
#                           'variants.map(v => va.calldata.combined.AN).stats().max']))


def subset_tj(hc):
    genomes = False

    vqsr_vds = None

    if(genomes):
        sample_lists = {'sa.concordance': genomes_concordance_}
        projects = '"G68758","G87944","G77318","G29747","G77318","G77525","G29749","G87944","G87944","G31561","G77318","G89387","G87944","G84381","G29747","G89387","G31561","G29748","G29748","G26842"'
        projects_expr = '[%s].toSet.contains(sa.project_or_cohort)' % projects
        path = full_genome_vds_path
        out = "genomes"
    else:
        sample_lists = {'sa.concordance': exomes_concordance_samples,
                        'sa.hapmap': 'gs://gnomad-exomes-raw/gnomad_exomes_hapmap_samples_ids.txt'}
        projects = '"C1629","FINNMETSEQ","FINNMETSEQ","C1975","C1476","UK10K","C1975","UK10K","C1476","C1975","UK10K","UK10K","C1753","C1975","C1622","FINNMETSEQ","C1629","C1975","C1629","C1975"'
        projects_expr = '[%s].toSet.contains(sa.meta.pid) || sa.hapmap' % projects
        path = full_exome_vds_path
        out = "exomes"
        vqsr_vds = hc.read(vqsr_vds_path)

    vds = hc.read(path)

    for root,sample_list in sample_lists.iteritems():
        vds = vds.annotate_samples_list(sample_list, root)

    filter_samples_expr = " || ".join(sample_lists.keys())

    if not genomes:
        annotations = ['culprit', 'POSITIVE_TRAIN_SITE', 'NEGATIVE_TRAIN_SITE', 'VQSLOD']
        vds = vds.annotate_variants_vds(vqsr_vds, code=', '.join(['va.info.%s = vds.info.%s' % (a, a) for a in annotations]))

    vds = (vds.annotate_samples_table(exomes_meta_tsv_path, 'sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
               .annotate_samples_expr(['sa.meta.project_description = sa.meta.description'])  # Could be cleaner
     )

    vds.export_samples('gs://gnomad-genomes/subsets/tj/exomes.sample_meta.txt.bgz', 'sa.meta.*')

    # vds = (
    #     vds
    #         .filter_samples_expr(filter_samples_expr)
    #         .filter_samples_expr(projects_expr)
    #         .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
    #         .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
    #         .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
    # )
    # print(vds.query_samples(['samples.count()']))
    # vds.export_vcf("gs://gnomad-genomes/subsets/tj/duplicate_samples.%s.vcf.bgz" % out)

#subset_tj(hc)

# paths = ["gs://gnomad/gnom.ad.vds","gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds"]
#
# for path in paths:
#     print("Running: %s\n\n" % path)
#     vds = (
#         hc.read(path,sites_only=True)
#             .annotate_variants_expr('va.v = v')
#             .split_multi()
#
#     )
#     print(vds.query_variants(['variants.map(v => v.start - va.v.start).counter()']))
#

#vds = hc.read(gnomad_path, sites_only=True)
# vds = hc.import_vcf("gs://gnomad-lfran/tmp/1var.vcf")
# #vds = vds.filter_variants_intervals(Interval.parse('1:1000000-1500000'))
# #vds = vds.filter_variants_intervals(Interval.parse('1:1087812-1087814'))
# #vds = vds.filter_variants_expr('v.altAlleles.forall(a => a.alt != "*")')
# vds = vds.vep(config=vep_config, csq=False, root='va.info.CSQ', force=True)
# #vds = vds.vep(config=vep_config, csq=False, root='va.info.CSQ', force=True)
# print(vds.variant_schema)
# print(vds.query_variants('variants.map(v => str(va.info.CSQ)).collect()'))
#
# vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)
# #vds = vds.vep(config=vep_config, csq=False, root='va.info.CSQ', force=True)
# print(vds.variant_schema)
# print(vds.query_variants('variants.map(v => str(va.info.CSQ)).collect()'))
# #vds.export_vcf("gs://gnomad-lfran/tmp/test-vep-attr.vcf")





# vds = hc.read(gnomad_path)
# vds = filter_intervals(vds,['2:21229160-21229160','2:21280079-21280079'])
# vds = vds.annotate_samples_expr('sa.keepMe = gs.filter(g => v.start == 21280079 && g.isCalledNonRef).count() > 0')
# vds = vds.filter_samples_expr('sa.keepMe')
# vds.export_vcf('gs://gnomad-lfran/tmp/APOB_Lof.vcf')

# y_intervals = '/tmp/chrY.txt'
# with open(y_intervals, 'w') as f:
#     f.write('Y:1-1000000000')
#
#
#
#
# all = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.X.bad-multi.vds")
#
# extra_fields = {
#     "va.info.BaseQRankSum" : []
# }
#
# multi = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.X.vds")
# multi = copy_schema_attributes(all, multi)
# multi = multi.annotate_variants_expr('va.pass = va.filters.isEmpty')
#
# #print("\n----- ALL -----\n")
# #print_schema_attributes(all)
# #print("\n----- Multi -----\n")
# #print_schema_attributes(multi)
#
# res = all.annotate_variants_vds(multi,'va = if(isMissing(vds)) va else vds')
# res.write("gs://gnomad-genomes/sites/internal/gnomad.genomes.sites.X.vds")


#
#     print("All filters")
#     vds = hc.read(p)
#     (
#         vds
#             .split_multi()
#             .variants_keytable().aggregate_by_key(
#             key_condition='type = if(v.altAllele.isSNP) "snv" else if(v.altAllele.isIndel) "indel" else "other",'
#                           'filtered = !va.info.AS_FilterStatus[va.aIndex-1].isEmpty || !va.filters.isEmpty',
#             agg_condition='n = va.count()')
#             .to_dataframe()
#             .show()
#     )


#print(vds.query_variants("variants.map(v => v.altAllele.isSNP).counter()"))

#run_sanity_checks(hc.read("gs://gnomad-exomes/sites/gnomad.exomes.sites.autosomes.vds"),['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS'])
# vds = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.X.vds")
# print(vds.query_variants(['variants.count()',
#                           'variants.filter(v => v.inXNonPar).count()']))
# print(vds.query_variants(['variants.map(v => range(v.nAltAlleles)'
#                     '.map(i => (va.info.Hom_AFR[i] + va.info.Hom_AMR[i] + va.info.Hom_ASJ[i] + va.info.Hom_EAS[i] + '
#                     'va.info.Hom_FIN[i] + va.info.Hom_NFE[i] + va.info.Hom_OTH[i] - va.info.Hom[i]).abs).max).stats()',
#                     'variants.map(v => range(v.nAltAlleles)'
#                     '.map(i => (va.info.AC_AFR[i] + va.info.AC_AMR[i] + va.info.AC_ASJ[i] + va.info.AC_EAS[i] + '
#                     'va.info.AC_FIN[i] + va.info.AC_NFE[i] + va.info.AC_OTH[i] - va.info.AC[i]).abs).max).stats()',
#                     'variants.map(v => (va.info.AN_AFR + va.info.AN_AMR + va.info.AN_ASJ + va.info.AN_EAS '
#                     '+ va.info.AN_FIN + va.info.AN_NFE + va.info.AN_OTH - va.info.AN).abs).stats()'
#                    ]))

#print(hc.read(gnomad_path).query_variants('variants.map(v => v.contig).counter()'))
#print hc.read('gs://gnomad-lfran/tmp/gnomad.sites.tmp.2017-02-13_21-40.vds').query_variants(['variants.filter(v => isDefined(va.info.%s)).count()' % x for x in ('DS', 'END', 'MQ0', 'MQ', 'RAW_MQ')])

#vds = vds.filter_variants_intervals('file://' + y_intervals)
#vds = vds.filter_variants_intervals('gs://gnomad-lfran/tmp/1gene.intervals')
#print(vds.query_variants('variants.map(v => isDefined(va.info.DS)).counter()'))
# s =hc.read("gs://gnomad-lfran/tmp/gnomad.sites.tmp.2017-02-13_21-40.vds").variant_schema
# print(s)
# rf = [x for x in s.fields if x.name == "info"][0]
# print(rf)
#for f in rf.fields:




#pprint(hc.read(gnomad_path, sites_only=True).query_variants('variants.filter(v => pcoin(0.001)).map(v => v.start).stats()'))
# vds = hc.read(vds_path, sites_only=True)
# vds = vds.filter_variants_expr('pcoin(0.1)')
# vds.export_vcf('file:///test/test.before.vcf')
# vds = hc.import_vcf('/test.before.new.vcf')
# vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)
# print(vds.variant_schema)
# vds.export_vcf('file:///test/test.after.vcf')
#vds.write("gs://gnomad-lfran/tmp/vep.test.vds")

# vds = hc.import_vcf('/test.before.new.vcf')
# vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)
# vds.export_vcf('file:///test/test.after.vcf')

#vds = hc.read("gs://gnomad-lfran/tmp/vep.test.vds")
#pprint(vds.query_variants('variants.filter(v => pcoin(0.1)).map(v => va.info.CSQ).collect()'))

    #.export_vcf("gs://gnomad-lfran/tmp/vep.test.vcf")
#pprint(vds.query_variants('variants .filter(v => v.altAlleles.exists(a => a.alt.length >1032)).collect()'))
#pprint(vds.query_variants('variants.map(v => v.altAlleles.map(a => a.alt.length).max).stats().max'))

#mvs = vds.query_variants('variants.filter(v => va.mendel.length > 0).collect()')[0]

#print(mvs[:3])


#label_criteria = ['va.transmitted_singleton && v.altAllele.isIndel', 'va.mills']
#pprint(dict(zip(label_criteria, vds.query_variants(['variants.filter(x => %s).count()' % x for x in label_criteria]))))


#
# hc.read("gs://gnomad/gnomad.raw_hardcalls.tmp.vds").print_schema()
#
# pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
#
# rf_features = ['va.variantType',
#                 'va.info.QD',
#                 'va.info.MQ',
#                 'va.info.MQRankSum',
#                 'va.info.FS',
#                 'va.info.SOR',
#                 'va.info.InbreedingCoeff',
#                 'va.info.ReadPosRankSum',
#                 'va.stats.raw.nrq_median',
#                 'va.stats.raw.ab_median',
#                 'va.stats.raw.dp_median',
#                 'va.stats.raw.gq_median']
#
# rf = (hc.read(rf_path,sites_only=True)
#         .annotate_global_expr_by_variant('global.variantsByType = index(variants.map(v => va.variantType).counter(),key)')
#     )
# rf.show_globals()
# (
#     rf.filter_variants_expr('pcoin([1.0,1000000 / global.variantsByType[va.variantType].count].min)', keep=True)
#     .export_variants('gs://gnomad-lfran/gnomad.features_for_median.txt',
#                    condition='v=v,'
#                              'variantType=va.variantType,'
#                              'mqranksum=va.info.MQRankSum,'
#                              'readposranksum=va.info.ReadPosRankSum,'
#                              'ab_median=va.stats.raw.ab_median')
# )
#

# global_expr = ['global.%s = variants.filter(x => isMissing(%s)).count()' % (a,a) for a in rf_features ]
#
# rf_kt =rf.variants_keytable()
#
# print(rf_kt.schema())
#
# rf_kt = rf_kt.flatten().select(['va.info.MQRankSum'])
#
# print(rf_kt.schema())
#
# rf_df = rf_kt.to_dataframe()
# rf_df.printSchema()
# rf_df.show()
#
# MQmedian = rf_df.approxQuantile("`va.info.MQRankSum`", [0.5], 0.05)
#
# print(MQmedian)


# (
#     hc.read(rf_path,sites_only=True)
#     .annotate_global_expr_by_variant(global_expr)
#     .show_globals()
#
# )


# (hc.read("gs://gnomad/gnomad.raw_hardcalls.vds",sites_only=True)
#  .filter_variants_expr('v.altAlleles.exists(a => a.isSNP) && pcoin(0.92)',
#                        keep=False)
#  .filter_variants_expr('v.altAlleles.exists(a => !a.isSNP) && pcoin(0.6)',
#                        keep=False)
#  .export_variants(output= "gs://gnomad-lfran/tmp/variantTypeCheck.txt.bgz",
#                     condition='v = v, alts = v.altAlleles.map(a => a.alt).mkString(","), variantType = va.variantType')
#  )

#intervals = "gs://gnomad-lfran/tmp/maryam_variants.txt"

#vds = set_vcf_filters(hc, "gs://gnomad-lfran/tmp/maryam_variants.txt", rf_path, rf_ann_root, rf_snv_cutoff, rf_indel_cutoff, filters = {}, filters_to_keep = [], tmp_path = '/tmp')

# x=  (
#     hc.read(vds_path)
#         .annotate_global_expr_by_sample('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
#         .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=pyhail.TextTableConfig(impute=True))
#         .annotate_samples_expr(['sa.meta.population = sa.meta.predicted_pop',
#                                 'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
#         .filter_samples_expr('!isMissing(sa.meta.predicted_pop)')  # Could be cleaner
#         .filter_variants_intervals('gs://gnomad-lfran/tmp/test.interval')
#         .histograms('va.info')
#         .export_samples(output='file:///tmp/out_sampes.txt', condition='s.id,sa.meta.project_description')
#         .export_variants(output='file:///tmp/out.txt',
#                          condition='gq_hist_all = va.info.GQ_HIST_ALL,'
#                                    'gq_hist_alt = va.info.GQ_HIST_ALT,'
#                                    'dp_hist_all = va.info.DP_HIST_ALL,'
#                                    'dp_hist_alt = va.info.DP_HIST_ALT,'
#                                    'ab_hist_all = va.info.AB_HIST_ALL,'
#                                    'ab_hist_alt = va.info.AB_HIST_ALT'
#                          )
# )
# print(x.print_schema())

# res = annotate_non_split_from_split(hc, non_split_vds_path=vds_path,
#                               split_vds=hc.read(rf_path),
#                               annotations=['va.RF1'],
#                               annotation_exp_out_path=tmp_RF_ann_out)
#
#
#
# res.write(tmp_vds)
# rf = hc.read(rf_path)
# (
#     hc.read(tmp_vds)
#     .annotate_variants_expr('va.alts = v.altAlleles.map(a => a.alt)')
#     .split_multi()
#     .annotate_variants_vds(rf,code='va.rf_ann = vds.RF1')
#     .annotate_variants_expr('va.issame = va.rf_ann == va.RF1[va.aIndex-1]')
#     .export_variants(tmp_RF_ann_exp, 'v=v,va.alts=va.alts, ai1 = va.aIndex, va.rf_ann.prob = va.rf_ann.probability["TP"], va.RF1 =  va.RF1,va.issame=va.issame')
# )

#.annotate_variants_expr('va.rf_match = va.rf_ann.probability["TP"] == va.RF1[va.aIndex - 1].probability["TP"]')
#va.RF1.prob = va.RF1[va.aIndex - 1].probability["TP"],

send_message(channels='@laurent', message='Test is done processing!')