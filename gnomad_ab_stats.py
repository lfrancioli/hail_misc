from hail import HailContext, Interval
from resources import get_gnomad_data
import re

def get_ab_stats(vds, outprefix):
    for interval, filter_males in [("1-22",False),("X",True)]:
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

hc = HailContext(log='/abstats.log')
vds = get_gnomad_data(hc, "genomes", release_samples=True)
get_ab_stats(vds, "gs://gnomad-lfran/tmp/gnomad.abstats")