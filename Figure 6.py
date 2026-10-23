from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M6'
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
INPUT_FILES = {'cep-np/07_validation/gse56081/gse56081_sample_scores.tsv': 'Inputs/0c5a0e0581 GSE56081 sample scores.tsv', 'cep-np/07_validation/gse70362/gse70362_np_donor_scores.tsv': 'Inputs/397e7068bd GSE70362 np donor scores.tsv', 'cep-np/07_validation/gse244889/gse244889_bulk_sample_scores.tsv': 'Inputs/921ab48c48 GSE244889 bulk sample scores.tsv', 'cep-np/07_validation/gse186542/sample_scores.tsv': 'Inputs/f51384c1fc sample scores.tsv'}
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
import math
import shutil
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
ROOT = (HERE / 'Source keys/cep-np/scripts/48_build_final_manuscript_figures.py').resolve().parents[1]
OUT = ROOT / 'manuscript_package_v1/figures/main'
BLUE = '#275D8C'
ORANGE = '#C85A0A'
RED = '#B4443E'

def style():
    sns.set_theme(style='whitegrid', context='paper')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'figure.dpi': 150, 'savefig.dpi': 400, 'axes.spines.top': False, 'axes.spines.right': False})

def label_panel(ax, label):
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

def save(fig, name, output_dir=OUT):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f'{name}.png', bbox_inches='tight', facecolor='white')
    fig.savefig(output_dir / f'{name}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)

def hedges_g_ci(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n1 = len(a)
    n0 = len(b)
    df = n1 + n0 - 2
    pooled = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / df)
    d = (a.mean() - b.mean()) / pooled
    j = 1 - 3 / (4 * df - 1)
    g = j * d
    var = (n1 + n0) / (n1 * n0) + g * g / (2 * df)
    se = math.sqrt(var)
    return (g, g - 1.96 * se, g + 1.96 * se)

def figure4(output_dir=OUT, output_name='Figure4_bulk_heterogeneity_across_four_cohorts'):
    np.random.seed(20261020)
    cohorts = []
    samples = []
    d = pd.read_csv(ROOT / '07_validation/GSE56081/GSE56081_sample_scores.tsv', sep='\t')
    a = d.loc[d.condition.eq('degenerated'), 'transition_z']
    b = d.loc[d.condition.eq('control'), 'transition_z']
    cohorts.append(('GSE56081', *hedges_g_ci(a, b), 0.00794))
    samples.append(pd.DataFrame({'cohort': 'GSE56081', 'group': d.condition.replace({'control': 'reference', 'degenerated': 'disease'}), 'score': d.transition_z}))
    d = pd.read_csv(ROOT / '07_validation/GSE70362/GSE70362_np_donor_scores.tsv', sep='\t')
    d = d.loc[d.grade_group.isin(['I-II', 'IV-V'])]
    a = d.loc[d.grade_group.eq('IV-V'), 'transition_z']
    b = d.loc[d.grade_group.eq('I-II'), 'transition_z']
    cohorts.append(('GSE70362', *hedges_g_ci(a, b), 0.1097))
    samples.append(pd.DataFrame({'cohort': 'GSE70362', 'group': d.grade_group.replace({'I-II': 'reference', 'IV-V': 'disease'}), 'score': d.transition_z}))
    d = pd.read_csv(ROOT / '07_validation/GSE244889/GSE244889_bulk_sample_scores.tsv', sep='\t')
    a = d.loc[d.group.eq('SDD'), 'transition_z']
    b = d.loc[d.group.eq('MDD'), 'transition_z']
    cohorts.append(('GSE244889 bulk', *hedges_g_ci(a, b), 0.7))
    samples.append(pd.DataFrame({'cohort': 'GSE244889 bulk', 'group': d.group.replace({'MDD': 'reference', 'SDD': 'disease'}), 'score': d.transition_z}))
    d = pd.read_csv(ROOT / '07_validation/GSE186542/sample_scores.tsv', sep='\t')
    a = d.loc[d.group.eq('advanced'), 'transition_z']
    b = d.loc[d.group.eq('early'), 'transition_z']
    cohorts.append(('GSE186542', *hedges_g_ci(a, b), 1.0))
    samples.append(pd.DataFrame({'cohort': 'GSE186542', 'group': d.group.replace({'early': 'reference', 'advanced': 'disease'}), 'score': d.transition_z}))
    eff = pd.DataFrame(cohorts, columns=['cohort', 'g', 'low', 'high', 'p'])
    sam = pd.concat(samples, ignore_index=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={'width_ratios': [1, 1.35]})
    ax = axes[0]
    label_panel(ax, 'A')
    y = np.arange(len(eff))[::-1]
    for i, r in enumerate(eff.itertuples()):
        col = BLUE if r.g > 0 else ORANGE
        ax.plot([r.low, r.high], [y[i], y[i]], color=col, lw=1.6)
        ax.plot(r.g, y[i], 'o', color=col, ms=7)
        if i == 0:
            ax.annotate(f'P={r.p:.3g}', (r.g, y[i]), xytext=(0, -12), textcoords='offset points', ha='center', va='top', fontsize=7)
        else:
            ax.annotate(f'P={r.p:.3g}', (r.g, y[i]), xytext=(0, 7), textcoords='offset points', ha='center', va='bottom', fontsize=7)
    ax.axvline(0, color='black', lw=0.9)
    ax.set_yticks(y, eff.cohort)
    ax.set_xlabel('Hedges g disease minus reference')
    ax.set_title('Frozen bulk effect sizes')
    ax = axes[1]
    label_panel(ax, 'B')
    order = eff.cohort.tolist()
    sns.boxplot(data=sam, x='cohort', y='score', hue='group', order=order, palette={'reference': '#9E9E9E', 'disease': ORANGE}, showfliers=False, ax=ax)
    sns.stripplot(data=sam, x='cohort', y='score', hue='group', order=order, palette={'reference': '#5B5B5B', 'disease': RED}, dodge=True, size=4, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title='', frameon=False)
    ax.axhline(0, color='black', lw=0.7)
    ax.set_xlabel('')
    ax.set_ylabel('frozen 25-gene transition score')
    ax.tick_params(axis='x', rotation=22)
    ax.set_title('Within-cohort sample distributions')
    fig.suptitle('The single-cell-defined state does not yield a transportable bulk degeneration marker', fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    save(fig, output_name, output_dir)

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def writepdf(page, path):
    writer = PdfWriter()
    writer.add_page(page)
    with path.open('wb') as f:
        writer.write(f)

def crop(source, destination, top=27):
    page = PdfReader(source).pages[0]
    page.mediabox.top = float(page.mediabox.top) - top
    page.cropbox.top = page.mediabox.top
    writepdf(page, destination)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Figure 6 raw.pdf'
    output = OUTPUT_DIR / 'Figure 6.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        style()
        fixed = ORIGINALREAD(HERE / 'Inputs/Bulk effect display.tsv', sep='\t')
        effects = iter(fixed[['g', 'low', 'high']].itertuples(index=False, name=None))
        hedges_g_ci = lambda a, b: next(effects)

        def local_save(fig, name, output_dir=None):
            rawsave(fig, raw, bbox_inches='tight', facecolor='white')
            plt.close(fig)
        save = local_save
        figure4()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    crop(raw, output)
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
