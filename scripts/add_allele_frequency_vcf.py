#!/usr/bin/env python3

import pysam
import argparse
import sys

def make_header(invcf):
    invcf.header.formats.add("ALT_AF", "A", "Float", "Alternate allele frequencies")
    invcf.header.formats.add("SUM_ALT", "1", "Integer", "Sum of observations from all alternate alleles")
    invcf.header.formats.add("SUM_ALT_AF", "1", "Float", "Sum of alternate allele frequencies")
    if 'AD' not in invcf.header.formats.keys():
        invcf.header.formats.add("AD", "2", "Integer", "Number of observation of REF and ALT alleles")
    
def add_tags(variant, sample, source_format):
    if source_format == 'delly':
        dr = variant.samples[sample]['DR']
        dv = variant.samples[sample]['DV']
        rr = variant.samples[sample]['RR']
        rv = variant.samples[sample]['RV']

        # If the variant is "IMPRECISE" we use the variant pairs to calculate
        # alt allele freq; if PRECISE we use the junctions. See code in
        # filter.h from delly:
        #     if (!precise) rVar = (float) dv[i] / (float) (dr[i] + dv[i]);
        #     else rVar = (float) rv[i] / (float) (rr[i] + rv[i]);
        assert type(variant.info['PRECISE']) == bool
        if variant.info['PRECISE']:
            ref = rr
            alt = rv
        else:
            ref = dr
            alt = dv
        variant.samples[sample]['AD'] = (ref, alt)
    elif source_format == 'freebayes':
        pass
    else:
        raise Exception('Invalid source format: "%s"' % source_format)
    ad = variant.samples[sample]['AD']
    if len(ad) == 0 or ad[0] is None:
        dp = 0
    else:
        dp = sum(ad)
    if dp > 0:
        alt_af = [x/float(dp) for x in ad[1:]]
        sumalt = sum(ad) - ad[0]
        variant.samples[sample]['ALT_AF'] = alt_af
        variant.samples[sample]['SUM_ALT'] = sumalt
        variant.samples[sample]['SUM_ALT_AF'] = sum(alt_af)


parser = argparse.ArgumentParser(description= 'Add tags to VCF')
parser.add_argument('invcf', default= '-', nargs= '?', help= 'Input VCF [%(default)s]')
parser.add_argument('--format', '-f', default= 'freebayes', choices= ['freebayes', 'delly'], help= 'Format (source) of VCF input [%(default)s]')
parser.add_argument('--version', action='version', version='%(prog)s 0.3.0')
args = parser.parse_args()


invcf = pysam.VariantFile(args.invcf, "r")
make_header(invcf)

outvcf = pysam.VariantFile('-', 'w', header= invcf.header)

for variant in invcf:
    try:
        for sample in variant.samples:
            add_tags(variant, sample, source_format= args.format)
    except:
        sys.stderr.write('Error in\n' + str(variant) + '\n')
        raise 
    outvcf.write(variant)
invcf.close()
