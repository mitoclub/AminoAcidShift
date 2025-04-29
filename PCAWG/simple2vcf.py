import sys
import tqdm
import pandas as pd

simple_file = sys.argv[1]
outpath = sys.argv[2]

# o = open(outpath, 'w')
# o.write('##fileformat=VCFv4.2\n')
# o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

print('##fileformat=VCFv4.2')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

data = []
i = 0
with open(simple_file) as fin:
    for line in tqdm.tqdm(fin):
        i += 1
        row = line.strip().split('\t')
        try:
            disease, sid, cohort, assembly, vartype, chr, start, end, ref, var, annotations = row
        except ValueError:
            disease, sid, cohort, assembly, vartype, chr, start, end, ref, var = row
        # except:
        #     print(line.strip().split('\t'))
        #     break
        
        data.append([chr, int(start), '.', ref, var, '.', '.', '.', '.'])
        # o.write(
        #     '\t'.join([chr, start, '.', ref, var, '.', '.', '.', '.']) + '\n'
        # )
        # if i > 100000:
        #     break

# o.close()

chr_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
chr_index = range(24)
chr_index_map = dict(zip(chr_order, chr_index))
print(chr_index_map)

df = pd.DataFrame(data)
df = df[(df[0] != 'MT') & (df[0] != 'M')]
df['chr_index'] = df[0].map(chr_index_map, )

assert df['chr_index'].isna().sum() == 0

# TODO subset = [0,1,4,5]
df = df.drop_duplicates(subset=[0,1,4]).sort_values(['chr_index', 1])
del df['chr_index']
df.to_csv(outpath, header=False, index=False, sep='\t')
