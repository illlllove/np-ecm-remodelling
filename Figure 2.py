from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M2'
PNG_DPI = 300
WORK_DIR = OUTPUT_DIR / 'Temporary files' / FIGURE_ID
WORK_DIR.mkdir(parents=True, exist_ok=True)
sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', str(OUTPUT_DIR / 'Matplotlib cache'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(OUTPUT_DIR / 'Numba cache'))
dependency_path = HERE / 'Environment/Python/Lib/site-packages'
if dependency_path.is_dir():
    sys.path.insert(0, str(dependency_path))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
ORIGINALREAD = pd.read_csv
ORIGINALSAVE = Figure.savefig
INPUT_FILES = {'cep-np/manuscript_package_bmc_genomics_v3/editorial/presentation inputs/cluster markers.tsv.gz': 'Inputs/0613a079f0 cluster markers.tsv.gz', 'cep-np/manuscript_package_bmc_genomics_v3/editorial/presentation inputs/codetection summary.tsv': 'Inputs/d03832b02b codetection summary.tsv', 'results/gse165722_independent_validation/04_phenotype/identity_fixed_panel_descriptive_display.tsv': 'Inputs/66f4d311dd identity fixed panel descriptive display.tsv'}
MANIFEST = {'canonicalRoot': str(HERE / 'Source keys/cep-np/manuscript package BMC Genomics v6')}

def norm(path):
    return str(Path(path)).replace('\\', '/').lower()

def cached_read(path, *args, **kwargs):
    if hasattr(path, 'read'):
        return ORIGINALREAD(path, *args, **kwargs)
    key = norm(path).split('/source keys/', 1)[-1]
    cached = INPUT_FILES.get(key)
    if cached is None:
        raise FileNotFoundError('No cached plot input: ' + str(path))
    return ORIGINALREAD(HERE / cached, *args, **kwargs)
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from PIL import Image
import seaborn as sns
PKG = (HERE / 'Source keys/cep-np/manuscript_package_BMC_Genomics_v4/editorial/Render integrated figures.py').resolve().parents[1]
ROOT = PKG.parent
BASE = ROOT / 'manuscript_package_BMC_Genomics_v3'
HOMEBASE = ROOT.parent
G165 = HOMEBASE / 'results/GSE165722_independent_validation'
MAIN = PKG / 'figures/main'
INPUT = BASE / 'editorial/presentation inputs'
BLUE = '#377DA8'
ORANGE = '#C77030'
SOURCES = {}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read(path, **kwargs):
    path = Path(path)
    SOURCES[str(path)] = sha(path)
    return pd.read_csv(path, **kwargs)

def label(ax, letter, x=-0.13, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=15, fontweight='bold', ha='left', va='bottom')

def finish(fig, name, folder=MAIN):
    for ext in ['pdf', 'png']:
        fig.savefig(folder / f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def common_style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10.5, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

def figure2():
    fig = plt.figure(figsize=(14.6, 10.8))
    grid = fig.add_gridspec(2, 2, left=0.095, right=0.98, bottom=0.105, top=0.94, wspace=0.4, hspace=0.56, height_ratios=[1.05, 1])
    axes = [fig.add_subplot(grid[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    for ax, letter in zip(axes, 'ABCD'):
        label(ax, letter)
    markers = read(INPUT / 'cluster markers.tsv.gz', sep='\t')
    genes = ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2', 'COL6A1', 'COL6A2', 'SPARC', 'POSTN', 'PCOLCE', 'MMP2', 'MMP14', 'CTHRC1', 'FAP', 'AEBP1']
    q = markers[markers.group.astype(str).isin(['5', '11']) & markers.names.isin(genes)]
    matrix = q.pivot(index='names', columns='group', values='logfoldchanges').reindex(genes)
    assert matrix.shape == (15, 2) and (not matrix.isna().any().any())
    matrix.columns = ['State 5', 'State 11']
    sns.heatmap(matrix, cmap='YlOrRd', linewidths=0.35, linecolor='white', ax=axes[0], cbar_kws={'label': 'Cluster log fold change', 'pad': 0.02})
    axes[0].set(xlabel='', ylabel='', title='GSE233666 fibrillar matrix and processing programme')
    axes[0].set_yticklabels(genes, rotation=0, fontstyle='italic', fontsize=8)
    d = read(G165 / '04_phenotype/identity_fixed_panel_descriptive_display.tsv', sep='\t')
    pops = ['Cluster 15', 'Cluster 16', 'Cluster 18', 'Other NP-related']
    genes2 = ['KRT8', 'KRT18', 'ACAN', 'COL2A1', 'SOX9', 'COMP', 'CILP', 'PRG4', 'COL1A1', 'COL3A1', 'PCOLCE', 'MMP2', 'PECAM1', 'RGS5', 'PDGFRB']
    q = d[d.population.isin(pops) & d.gene.isin(genes2)].set_index(['population', 'gene'])
    rows = []
    for yi, pop in enumerate(pops):
        for xi, gene in enumerate(genes2):
            row = q.loc[pop, gene]
            rows.append((xi, yi, float(row.fraction), float(row.mean_logCP10K)))
    arr = np.array(rows)
    ax = axes[1]
    dots = ax.scatter(arr[:, 0], arr[:, 1], s=arr[:, 2] * 180, c=arr[:, 3], cmap='viridis', vmin=0, vmax=5, edgecolors='none')
    ax.set(xlim=(-0.7, 14.7), ylim=(3.55, -0.8), yticks=range(4), yticklabels=pops, title='GSE165722 conditional supplied-count audit')
    ax.set_xticks(range(15), labels=genes2, rotation=65, ha='right', fontstyle='italic', fontsize=8)
    ax.axvline(7.5, color='#C8CDD2', lw=0.9)
    ax.axvline(11.5, color='#C8CDD2', lw=0.9)
    for center, name in [(3.5, 'NP/cartilage'), (9.5, 'ECM'), (13, 'Vascular/mural')]:
        ax.text(center, -0.47, name, fontsize=7.5, ha='center')
    cb = fig.colorbar(dots, ax=ax, fraction=0.035, pad=0.025)
    cb.set_label('Mean log1p CP10K', fontsize=8)
    handles = [ax.scatter([], [], s=x * 180, c='#536B72', edgecolors='none', label=f'{int(x * 100)}%') for x in [0.1, 0.5, 1]]
    ax.legend(handles=handles, title='Detected cells', frameon=False, ncol=3, loc='upper left', bbox_to_anchor=(-0.02, -0.29), fontsize=7.5, title_fontsize=8, handletextpad=0.3, columnspacing=1.2)
    data = read(INPUT / 'codetection summary.tsv', sep='\t')
    contexts = ['GSE244889', 'GSE233666', 'GSE222182']
    labels = ['GSE244889\n7 donors', 'GSE233666\n4 patients', 'GSE222182\n6 samples']
    data = data[data.population.eq('Frozen remodelling')].set_index(['context', 'event', 'threshold_umi'])
    for ax, events, colours, legendlabels, title in [(axes[2], [('COL1A1_cartilage_ge1', 1), ('COL1A1_cartilage_ge1', 2)], [BLUE, ORANGE], ['At least 1 UMI per gene', 'At least 2 UMIs per gene'], 'Fibrous and cartilage transcripts in the same cells'), (axes[3], [('FAP_COL1A1', 1), ('SFRP4_COL1A1', 1)], ['#5B8E73', '#9B6A9E'], ['$\\mathit{FAP}$ + $\\mathit{COL1A1}$', '$\\mathit{SFRP4}$ + $\\mathit{COL1A1}$'], 'Variable subtype-associated expression backgrounds')]:
        for k, (event, threshold) in enumerate(events):
            selected = data.loc[[(c, event, threshold) for c in contexts]]
            med = selected.median_fraction.to_numpy() * 100
            errors = np.vstack([med - selected.q1_fraction.to_numpy() * 100, selected.q3_fraction.to_numpy() * 100 - med])
            x = np.arange(3) + (k - 0.5) * 0.32
            ax.bar(x, med, width=0.29, color=colours[k], label=legendlabels[k], alpha=0.94)
            ax.errorbar(x, med, yerr=errors, fmt='none', ecolor='#333333', capsize=3, linewidth=1)
            for xx, yy, top in zip(x, med, selected.q3_fraction.to_numpy() * 100):
                ax.text(xx, top + 1.4, f'{yy:.1f}', ha='center', va='bottom', fontsize=8)
        ax.set_xticks(np.arange(3), labels=labels)
        ax.set(ylabel='Median positive state cells (%)', ylim=(0, 105), title=title)
        ax.legend(frameon=False, loc='upper left', fontsize=8)
        ax.grid(axis='y', alpha=0.16)
        ax.set_axisbelow(True)
    fig.text(0.095, 0.038, 'C/D: unchanged donor/sample medians and IQRs. Cartilage: $\\mathit{ACAN}$, $\\mathit{COL2A1}$ or $\\mathit{SOX9}$ with $\\mathit{COL1A1}$.', fontsize=8)
    fig.text(0.095, 0.018, 'B is descriptive independent-cohort support, not grade validation or evidence that all candidates form one stable NP subtype.', fontsize=8)
    finish(fig, 'Figure 2')

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Figure 2 raw.pdf'
    output = OUTPUT_DIR / 'Figure 2.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        common_style()
        read = cached_read
        figure2()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    output.write_bytes(raw.read_bytes())
    import pypdfium2 as pdfium
    document = pdfium.PdfDocument(str(output))
    page = document[0]
    bitmap = page.render(scale=PNG_DPI / 72)
    image = bitmap.to_pil()
    image.save(output.with_suffix('.png'))
    image.close()
    bitmap.close()
    page.close()
    document.close()
    return output
if __name__ == '__main__':
    make_figure()
