import pandas as pd
import snakemake as smake

lens = list()
shell_script = "zcat {fastq} | wc -l || true"

for dataset, fastq in snakemake.input.items():
    for size in smake.shell(shell_script.format(fastq=fastq), iterable=True):
        lens.append((dataset, 1/(int(size)/1e6)))

df = pd.DataFrame(lens, columns=['dataset', 'size'])
df.to_csv(snakemake.output.datasets_depth, sep='\t', index=None)
