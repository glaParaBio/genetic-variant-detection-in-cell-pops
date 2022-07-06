#!/usr/bin/env python3

import sys
import argparse
from collections import OrderedDict
from urllib.parse import unquote

parser = argparse.ArgumentParser(description= 'Extract gene_id, Name, and description from GFF')

parser.add_argument('gff', help= 'Input GFF file or stdin [%(default)s]', default= '-', nargs= '?')
parser.add_argument('--version', '-v', action= 'version', version= '%(prog)s 0.1.0')

args = parser.parse_args()

if args.gff == '-':
    gff = sys.stdin
else:
    gff = open(args.gff)

attribute_keys = ['Name', 'description']
header = ['gene_id', 'gene_name', 'description']
print('\t'.join(header))

table = OrderedDict()

for line in gff:
    if line.startswith('#'):
        continue
   
    attr_list = line.strip().split('\t')[8].split(';')
    attributes = {}
    for x in attr_list:
        k,v = x.split('=')
        attributes[k] = v
    
    if all([x not in attributes for x in attribute_keys]):
        continue
    
    gene_id = None
    for gene_id_key in ['gene_id', 'Parent', 'ID']:
        if gene_id_key in attributes:
            gene_id = attributes[gene_id_key]
            break
    assert gene_id is not None
    
    if gene_id not in table:
        table[gene_id] = {}
        for x in attribute_keys:
            table[gene_id][x] = ''

    for x in attribute_keys:
        if x not in attributes:
            continue
        
        if x == 'description':
            attributes[x] = unquote(attributes[x])

        if table[gene_id][x] == '':
            table[gene_id][x] = attributes[x]
        else:
            assert table[gene_id][x] == attributes[x]

for gene_id in table:
    outline = [gene_id]
    for x in attribute_keys:
        outline.append(table[gene_id][x])
    print('\t'.join(outline))

gff.close()
