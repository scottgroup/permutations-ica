import pandas as pd


def get_col_tuple(file):
    return tuple([col for col in file.split('/')[2:-1] if col not in ('stranded', 'unstranded')])


def parse_cufflinks(file):
    _df = pd.read_csv(file, sep=' ', header=None, index_col=0)
    return _df[1].to_dict()


def parse_featureCounts(file):
    _df = pd.read_csv(
        file, sep=' ',
        names=['len', 'count'], header=0, index_col=0
    )
    return _df['count'].to_dict()


def parse_htseq(file):
    _df = pd.read_csv(file, sep='\t', header=None, index_col=0)
    return _df[1].to_dict()


# Extracting count results from the different quantifiers
annot = snakemake.wildcards.annotation
if annot =='refseq':
    names = ['symbol', 'gene', 'HGNC_ID']
else:
    names = ['symbol', 'gene']

files = [file for file in snakemake.input if annot in file]
data_df = pd.read_csv(
    snakemake.input.gene_list,
    sep='\t', header=None, names=names
)

if annot == 'refseq':
    data_df.set_index('symbol', inplace=True)
else:
    data_df.set_index('gene', inplace=True)

for file in files:
    file_dict = ''
    if 'cufflinks' in file:
        file_dict = parse_cufflinks(file)
    elif 'featureCounts' in file:
        file_dict = parse_featureCounts(file)
    elif 'htseq' in file:
        file_dict = parse_htseq(file)

    if file_dict:
        col_tuple = get_col_tuple(file)
        data_df[col_tuple] = data_df.index.map(file_dict)

data_df.reset_index(inplace=True)
data_df.set_index(names, inplace=True)

print(data_df)
print(data_df.columns)

data_df.columns = pd.MultiIndex.from_tuples(
    data_df.columns,
    names=['quantifier','tissue', 'data_id', 'trimmer', 'aligner', 'annotation']
)

# Saving to file
data_df.to_csv(snakemake.output.out_NaN_file, sep='\t')

# Changing NaN to zero
data_df.fillna(value='F', inplace=True)

# Saving no NaN to file
data_df.to_csv(snakemake.output.out_file, sep='\t')
