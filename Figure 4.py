from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M4'
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
INPUT_FILES = {'cep-np/09_directiona_composition_bulk/02_differential_abundance/formal_abundance_effects.tsv': 'Inputs/2a5eb793db formal abundance effects.tsv', 'cep-np/07_validation/gse244889/formal_donor_effects.tsv': 'Inputs/5c6934a63c formal donor effects.tsv', 'cep-np/09_frozen_tables/wp1_wp5_external_validation_v1/gse244889_strict_marker_identity_sensitivity.tsv': 'Inputs/ab129d0717 GSE244889 strict marker identity sensitivity.tsv', 'cep-np/09_directiona_composition_bulk/03_state_pseudobulk/remodeling_state_score_effects.tsv': 'Inputs/c36d96a084 remodeling state score effects.tsv', '独立审查/表型与供体审查 20261021/gse244889 donor audit.csv': 'Inputs/977ed4272a GSE244889 donor audit.csv', '独立审查/表型与供体审查 20261021/gse244889 leave one donor audit.csv': 'Inputs/34d7e30b7d GSE244889 leave one donor audit.csv'}
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
HOMEBASE = ROOT.parent
AUDIT = HOMEBASE / '独立审查/表型与供体审查 20261021'
MAIN = PKG / 'figures/main'
BLUE = '#377DA8'
ORANGE = '#C77030'
GREEN = '#498F84'
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

def abundance_data():
    a = read(ROOT / '09_directionA_composition_bulk/02_differential_abundance/formal_abundance_effects.tsv', sep='\t')
    g = read(ROOT / '07_validation/GSE244889/formal_donor_effects.tsv', sep='\t')
    strict = read(ROOT / '09_frozen_tables/WP1_WP5_external_validation_v1/GSE244889_strict_marker_identity_sensitivity.tsv', sep='\t')
    sdiff = strict.loc[strict.condition.eq('SDD'), 'fraction'].mean() - strict.loc[strict.condition.eq('MDD'), 'fraction'].mean()
    rows = [[r.study_id, r.difference, r.bootstrap_ci_low, r.bootstrap_ci_high, r.exact_p] for r in a.itertuples()]
    r = g[g.analysis.eq('remodeling_fraction')].iloc[0]
    rows.append(['GSE244889', r.SDD_minus_MDD, r.ci_low, r.ci_high, r.exact_p])
    rows.append(['GSE244889 strict', sdiff, np.nan, np.nan, 1 / 35])
    return pd.DataFrame(rows, columns=['cohort', 'effect', 'low', 'high', 'p'])

def forest(ax, d, color=BLUE, percent=False):
    scale = 100 if percent else 1
    for i, r in enumerate(d.itertuples()):
        col = GREEN if 'strict' in r.cohort else color
        ax.scatter(r.effect * scale, i, s=30, c=col, zorder=3)
        if np.isfinite(r.low):
            ax.plot([r.low * scale, r.high * scale], [i, i], color=col, lw=1.5)
        ax.annotate(f'P={r.p:.3g}', (r.effect * scale, i), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=7)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_yticks(range(len(d)), labels=d.cohort)
    ax.set_ylim(len(d) - 0.5, -0.7)
    ax.grid(axis='x', alpha=0.16)
    ax.set_axisbelow(True)

def figure4():
    fig = plt.figure(figsize=(15.5, 10.8))
    grid = fig.add_gridspec(2, 2, left=0.11, right=0.975, bottom=0.09, top=0.94, wspace=0.44, hspace=0.48)
    sub = grid[0, 0].subgridspec(1, 2, width_ratios=[1.1, 1], wspace=0.48)
    aa = fig.add_subplot(sub[0])
    ab = fig.add_subplot(sub[1])
    b = fig.add_subplot(grid[0, 1])
    c = fig.add_subplot(grid[1, 0])
    d = fig.add_subplot(grid[1, 1])
    label(aa, 'A', x=-0.55)
    label(b, 'B')
    label(c, 'C')
    label(d, 'D')
    forest(aa, abundance_data(), percent=True)
    aa.set(title='Abundance effects', xlabel='Disease minus reference\n(percentage points)')
    donor = read(AUDIT / 'GSE244889 donor audit.csv')
    offsets = np.array([-0.15, -0.04, 0.1, 0.1, -0.08, 0, 0.08])
    donor_labels = {'Pa-17F': (-0.32, 2.6), 'Pab-30M': (0.12, 6.5), 'Pab-24M': (-0.02, 10.1), 'Pb-55F': (0.13, 14.8)}
    for index, r in enumerate(donor.itertuples()):
        x = (r.condition == 'SDD') + offsets[index]
        col = ORANGE if r.condition == 'SDD' else BLUE
        ab.scatter(x, r.remodeling_fraction * 100, c=col, s=30, zorder=3)
        if r.donor_id in donor_labels:
            ab.annotate(r.donor_id, (x, r.remodeling_fraction * 100), xytext=donor_labels[r.donor_id], textcoords='data', fontsize=6.8, arrowprops={'arrowstyle': '-', 'color': '#A8AFB5', 'lw': 0.5})
        else:
            ab.annotate(r.donor_id, (x, r.remodeling_fraction * 100), xytext=(3, 4), textcoords='offset points', fontsize=6.8)
    ab.set(xticks=[0, 1], xticklabels=['Mild\n4 donors', 'Severe\n3 donors'], xlim=(-0.45, 1.65), ylim=(0, 63), ylabel='State fraction (%)', title='Every GSE244889 donor')
    ab.grid(axis='y', alpha=0.16)
    ab.set_axisbelow(True)
    scores = read(ROOT / '09_directionA_composition_bulk/03_state_pseudobulk/remodeling_state_score_effects.tsv', sep='\t').query('threshold==30')
    vals = [[r.study_id, r.difference, r.bootstrap_ci_low, r.bootstrap_ci_high, r.exact_p] for r in scores.itertuples()]
    x = read(ROOT / '07_validation/GSE244889/formal_donor_effects.tsv', sep='\t').query("analysis=='within_state_score_n30'").iloc[0]
    vals.append(['GSE244889', x.SDD_minus_MDD, x.ci_low, x.ci_high, x.exact_p])
    forest(b, pd.DataFrame(vals, columns=['cohort', 'effect', 'low', 'high', 'p']), color=ORANGE)
    b.set(title='Additional within-state expression remains uncertain', xlabel='Disease minus reference frozen score')
    offset = {'Pa-17F': (-4, -14), 'Pab-24M': (4, 6), 'Pab-30M': (4, -13), 'Pb-55F': (-45, 6), 'Pc-41M': (5, 5), 'Pd-59F': (-45, 5), 'Pd-62F': (-45, 5)}
    for r in donor.itertuples():
        col = ORANGE if r.condition == 'SDD' else BLUE
        c.scatter(r.age, r.remodeling_fraction * 100, c=col, s=38, zorder=3)
        c.annotate(r.donor_id, (r.age, r.remodeling_fraction * 100), xytext=offset[r.donor_id], textcoords='offset points', fontsize=8)
    c.set(xlabel='Age (years)', ylabel='State fraction (%)', xlim=(12, 70), ylim=(0, 65), title='GSE244889 age context, not an adjusted disease effect')
    c.grid(alpha=0.16)
    c.legend(handles=[Line2D([], [], marker='o', linestyle='', color=BLUE, label='Mild: median age 27'), Line2D([], [], marker='o', linestyle='', color=ORANGE, label='Severe: median age 59')], frameon=False, loc='upper left', fontsize=8)
    lodo = read(AUDIT / 'GSE244889 leave one donor audit.csv')
    yy = np.arange(len(lodo))
    d.scatter(lodo.effect * 100, yy, c='#364554', s=30, zorder=3)
    d.axvline(26.3732794862, color='#777777', lw=1, ls='--')
    d.axvline(0, color='black', lw=0.8)
    for i, r in enumerate(lodo.itertuples()):
        d.text(38, i, f'P={r.exactP:.4f}', va='center', fontsize=8)
    d.set(yticks=yy, yticklabels=lodo.omitted, xlim=(0, 49), ylim=(len(lodo) - 0.5, -0.5), xlabel='Severe minus mild (percentage points)', title='LODO direction persists, significance does not')
    d.grid(axis='x', alpha=0.16)
    d.set_axisbelow(True)
    fig.text(0.11, 0.025, 'A/B retain frozen primary effects. C is descriptive. D exact P values were newly enumerated in the independent donor audit; no biological labels were changed.', fontsize=8)
    finish(fig, 'Figure 4')

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Figure 4 raw.pdf'
    output = OUTPUT_DIR / 'Figure 4.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        common_style()
        read = cached_read
        figure4()
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
