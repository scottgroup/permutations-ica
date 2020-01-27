import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

from plotting_utils import mm2inch, rcParams, colors
for k, v in rcParams.items():
    plt.rcParams[k] = v


def get_mult_align(exon, seq):
    """ Get multi alignment """
    hits = list()
    a = pairwise2.align.localms(exon, seq, 1, -0.75, -2, -0.25)[0]
    if isinstance(a[0], list):
        a = a[0]

    _exon, _seq, score, start, end = a
    exon_start, exon_end = 0, len(exon)
    seq_start, seq_end = 0, len(_seq)

    # If seq begins -
    while _seq[exon_start] == '-':
        exon_start += 1

    # If seq ends with -
    n_trailing = 0
    while _seq[seq_end-1] == '-':
        seq_end -= 1
        n_trailing += 1
    exon_end -= (len(_seq)- seq_end)
    _exon = _exon[exon_start:len(_exon)-n_trailing]

    if score > (exon_end-exon_start)/1.5:
        exon = exon[exon_start:exon_end]
        seq = _seq[exon_start:seq_end]

        for sub_exon in [sub for sub in _exon.split('-') if sub]:
            _ex_start = exon.find(sub_exon) + exon_start
            pos = _exon.find(sub_exon)
            mod_sub_exon = seq[pos:pos+len(sub_exon)]
            _seq_start = seq.find(mod_sub_exon)

            hits.append([_ex_start, _ex_start+len(sub_exon), _seq_start, _seq_start+len(sub_exon)])

    new_hits = list()
    gapped = 0
    for i, hit in enumerate(hits):
        _ex_start, _ex_end, _seq_start, _seq_end = hit
        local_seq = seq[_seq_start:_seq_end]
        if '-' in local_seq:
            sub_seqs = [sub for sub in local_seq.split('-') if sub]
            gapped += (len(local_seq) - len(''.join(sub_seqs)))
            running_len = 0
            for sub_seq in sub_seqs:
                pos = local_seq.find(sub_seq)
                new_hits.append([_ex_start+running_len, _ex_start+running_len+len(sub_seq), _seq_start+pos, _seq_start+pos+len(sub_seq)])
                running_len += len(sub_seq)
        else:
            new_hits.append(hit)

    # Correcting for gapped removed
    if new_hits:
        new_hits[-1][3] = new_hits[-1][3] - gapped
        new_hits[-1][1] = new_hits[-1][1] - gapped
    return(new_hits)


# Plotting!
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=mm2inch((178, 60)), gridspec_kw={'height_ratios':[4,1,4]}
)

# Read gene seq
for seq_record in SeqIO.parse(snakemake.input.gene_seq, "fasta"):
        gene_seq = seq_record
for seq_record in SeqIO.parse(snakemake.input.pseudo_seq, "fasta"):
        pseudo_seq = seq_record

# Loading GTF file
transcript_id = snakemake.wildcards.transcript_id
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4],
    names=[
        'seqname', 'source', 'feature', 'start',
        'end', 'score', 'strand', 'frame', 'attributes'
    ]
)
gtf['gene_id'] = gtf['attributes'].str.extract("gene_id \"(.*?)\"")
gtf['transcript_id'] = gtf['attributes'].str.extract("transcript_id \"(.*?)\"")
gtf_gene = gtf[gtf['gene_id'] == snakemake.wildcards.gene]
gtf_transcript = gtf_gene[(gtf_gene['transcript_id'] == transcript_id) & (gtf_gene['feature'] == 'exon')]
gene_line = gtf_gene[gtf_gene['feature'] == 'gene']
gene_start, gene_end = int(gene_line['start']), int(gene_line['end'])
gene_strand = str(gene_line['strand'].values[0])

gtf_pseudo = gtf[gtf['gene_id'] == snakemake.wildcards.pseudogene]
pseudo_line = gtf_pseudo[gtf_pseudo['feature'] == 'gene']
pseudo_strand = str(pseudo_line['strand'].values[0])

if gene_strand == '-':
    gene_seq = gene_seq.reverse_complement()
gene_seq = str(gene_seq.seq)

if pseudo_strand == '-':
    pseudo_seq = pseudo_seq.reverse_complement()
pseudo_seq = str(pseudo_seq.seq)

# Count files
gene = pd.read_csv(snakemake.input.gene, sep='\t', index_col=0, header=[0, 1, 2, 3, 4])
gene = gene.groupby(level=['aligner', 'dataset', 'tissue', 'trimmer'], axis=1).mean()
gene = gene.groupby(level='aligner', axis=1).mean()
aligners = list(gene.columns)

# Count files
pseudo = pd.read_csv(snakemake.input.pseudo, sep='\t', index_col=0, header=[0, 1, 2, 3, 4])
pseudo = pseudo.groupby(level=['aligner', 'dataset', 'tissue', 'trimmer'], axis=1).mean()
pseudo = pseudo.groupby(level='aligner', axis=1).mean()

# Parameters
extra = 0
padding = 50
n_exon = len(gtf_transcript)
len_vector = np.sum(gtf_transcript['end']+1-gtf_transcript['start']) + 2*extra*n_exon + padding*(n_exon-1)
vector = {aligner:np.zeros(len_vector) for aligner in aligners}

# Need to pad to counts
pseudo_counts = {aligner:np.zeros(len_vector) for aligner in aligners}
pad_pseudo = int((len_vector-len(pseudo))/4)
for aligner in aligners:
    pseudo_counts[aligner][pad_pseudo:(pad_pseudo+len(pseudo))] = pseudo[aligner]
# Drawing Gene structure
axes[1].add_artist(Rectangle(
    xy=(pad_pseudo, 0), width=len(pseudo), height=0.25, edgecolor=None
))


# Calculating exon distribution
x_pos = 0
exons = list()
for i, (idx, exon) in enumerate(gtf_transcript.iterrows()):
    vec_start = x_pos
    vec_end = x_pos+2*extra+(exon['end']+1-exon['start'])
    exon_start = exon['start']-gene_start-extra
    exon_end = exon['end']+1-gene_start+extra
    for aligner in aligners:
        _seq = gene[aligner].values[exon_start:exon_end]
        if gene_strand == '-':
            _seq = _seq[::-1]
        vector[aligner][vec_start:vec_end] = _seq
    exons.append((vec_start, vec_end))
    x_pos += 2*extra + padding + exon['end']-exon['start'] + 1

    # Calculating sequence alignment
    if gene_strand == '-':
        _start = exon_start
        exon_start = len(gene_seq) - exon_end
        exon_end = len(gene_seq) - _start
    _exon = gene_seq[exon_start:exon_end]
    hits = get_mult_align(_exon, str(pseudo_seq))

    # Draw ALL THE LINES
    top = 0.75
    bot = 0.25
    if len(hits) > 0:
        for _ex_start, _ex_end, _seq_start, _seq_end in hits:
            points = list()
            points.append((_ex_start+vec_start, top))
            points.append((_seq_start+pad_pseudo, bot))
            points.append((_seq_end+pad_pseudo, bot))
            points.append((_ex_end+vec_start, top))
            axes[1].add_artist(Polygon(
                xy=points, alpha=0.3, edgecolor=None
            ))

            # Adding vertical line in bottom
            axes[1].add_line(Line2D(
                [_seq_start+pad_pseudo, _seq_start+pad_pseudo], [0, 0.245], color='gray', linestyle='--', linewidth=0.5
            ))
            axes[1].add_line(Line2D(
                [_seq_end+pad_pseudo, _seq_end+pad_pseudo], [0, 0.245], color='gray', linestyle='--', linewidth=0.5
            ))

            # Add verical lines on the bedgraphs!
            max_start = np.max(pseudo.loc[_seq_start])
            axes[2].add_line(Line2D(
                [_seq_start+pad_pseudo, _seq_start+pad_pseudo], [0, max_start], color='gray', linestyle='--', linewidth=0.5
            ))


    # Drawing Gene structure
    axes[1].add_artist(Rectangle(
        xy=(vec_start, 0.75), width=vec_end-vec_start, height=0.25, edgecolor=None
    ))


# Drawing gene BedGraph
x_vec = np.arange(len(vector[aligners[0]]))
for aligner in aligners:
    for exon in exons:
        _first = exon[0]-1
        _second = exon[1]+2
        if _first < 0:
            _first = 0
        axes[0].plot(x_vec[_first:_second], vector[aligner][_first:_second], color=colors[aligner], linewidth=2, solid_capstyle='round')
    axes[2].plot(
        x_vec[pad_pseudo:len(pseudo)+pad_pseudo], pseudo_counts[aligner][pad_pseudo:len(pseudo)+pad_pseudo],
        color=colors[aligner], linewidth=2
    )

# Gene must have the same axis
xlim_gene = axes[0].get_xlim()
axes[1].set_xlim(xlim_gene)
axes[2].set_xlim(xlim_gene)
axes[0].axis('off')
axes[1].axis('off')
axes[2].axis('off')

ylim_flip = axes[2].get_ylim()
axes[2].set_ylim(ylim_flip[::-1])

plt.subplots_adjust(wspace=0, hspace=0, top=1, bottom=0, left=0.1, right=1)
plt.savefig(snakemake.output.plot)
