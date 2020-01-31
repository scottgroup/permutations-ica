import gzip

#  snakemake
from snakemake.shell import shell

chr_dict = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT"
}

gtf_file = shell(
    """
    zcat {snakemake.input}
    """,
    iterable=True
)

gtf = list()
for line in gtf_file:
    _line = line.split('\t')
    _chr = _line[0]
    if line[0] == '#':
        gtf.append(line)
    if _chr in chr_dict.keys():
        _line[0] = chr_dict[_chr]
        gtf.append('\t'.join(_line))

with gzip.open(snakemake.input[0] + '.temp', 'wb') as f:
    for line in gtf:
        f.write(line.encode() + '\n'.encode())

shell(
    "touch {snakemake.output}; "
    "rm {snakemake.input}; "
    "mv {snakemake.input}.temp {snakemake.input}"
)
