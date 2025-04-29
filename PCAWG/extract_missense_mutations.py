import pandas as pd


def read_annot(path):
    _cols = 'chr start ID ref var QUAL FILTER INFO FORMAT'.split()
    df = pd.read_csv(
        path, sep='\t', comment='#', header=None, names=_cols, 
        usecols=['chr', 'start', 'ref', 'var', 'INFO'],
        dtype={'chr': str},
    )
    df = df[df.INFO != '.']
    df = df[df.INFO.str.contains('missense_variant')]
    df['aa_change'] = df['INFO'].apply(lambda x: x.split('|')[15])
    df['cdn_change'] = df['INFO'].apply(lambda x: x.split('|')[16])
    df['strand'] = df['INFO'].apply(lambda x: x.split('|')[19])
    del df['INFO']
    return df

wgs38_annot = read_annot('./data/WGS.SNP.GRCh38.annotated.vcf')
wgs37_annot = read_annot('./data/WGS.SNP.GRCh37.annotated.vcf')
wes38_annot = read_annot('./data/WES.SNP.GRCh38.annotated.vcf')
wes37_annot = read_annot('./data/WES.SNP.GRCh37.annotated.vcf')

print('wgs38_annot')
print(wgs38_annot)

_names = 'disease sample cohort assembly snp/indel chr start end ref var annotations'.split()

wes38_mut = pd.read_csv('./data/WES.SNP.GRCh38.simple', sep='\t', header=None, names=_names, dtype={'chr': str},)\
    .assign(dataset='WES')
wes37_mut = pd.read_csv('./data/WES.SNP.GRCh37.simple', sep='\t', header=None, names=_names, dtype={'chr': str},)\
    .assign(dataset='WES')
wgs38_mut = pd.read_csv('./data/WGS.SNP.GRCh38.simple', sep='\t', header=None, names=_names, dtype={'chr': str},)\
    .assign(dataset='WGS')
wgs37_mut = pd.read_csv('./data/WGS.SNP.GRCh37.simple', sep='\t', header=None, names=_names, dtype={'chr': str},)\
    .assign(dataset='WGS')

wes38_mut_annot = wes38_mut.merge(wes38_annot)
wes37_mut_annot = wes37_mut.merge(wes37_annot)
wgs38_mut_annot = wgs38_mut.merge(wgs38_annot)
wgs37_mut_annot = wgs37_mut.merge(wgs37_annot)

print(len(wes38_annot), len(wes37_annot), len(wgs38_annot), len(wgs37_annot))

mut_total = pd.concat([wes37_mut_annot, wes38_mut_annot, 
                       wgs37_mut_annot, wgs38_mut_annot], ignore_index=True)

mut_total.to_csv('./data/missense_mutations.csv', index=False)

print(mut_total.groupby(['dataset', 'assembly']).chr.count().unstack())

print(mut_total)
