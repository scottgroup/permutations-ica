import pandas as pd
from snakemake import shell

resp = shell(
    "cat {snakemake.input.gtf} | "
    "grep -Po 'db_xref \"GeneID:\K.*?(?=\")' | "
    "sort | uniq -c"
    , iterable=True
)

data = pd.DataFrame([line.split() for line in resp], columns=['N', 'GeneID'])
data.to_csv(snakemake.output.GeneID_numbers, sep='\t', index=False)
