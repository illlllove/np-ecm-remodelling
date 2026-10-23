from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S22'
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
INPUT_FILES = {'cep-np/09_directiona_composition_bulk/04_decomposition/symmetric_coordinate_decomposition.tsv': 'Inputs/14f17ad232 symmetric coordinate decomposition.tsv', '独立审查/表型与供体审查 20261021/all donor proportions.csv': 'Inputs/d7af3ad9a1 All donor proportions.csv'}
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
SUPP = PKG / 'figures/supplementary'
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

def supplementary22():
    fig = plt.figure(figsize=(14.8, 9.8))
    grid = fig.add_gridspec(2, 3, left=0.095, right=0.97, bottom=0.11, top=0.94, wspace=0.5, hspace=0.55, height_ratios=[1, 1.5])
    ax = fig.add_subplot(grid[0, :])
    label(ax, 'A', x=-0.04)
    dec = read(ROOT / '09_directionA_composition_bulk/04_decomposition/symmetric_coordinate_decomposition.tsv', sep='\t')
    xx = np.arange(len(dec))
    w = 0.25
    ax.bar(xx - w / 2, dec.composition_component, w, color=BLUE, label='Composition')
    ax.bar(xx + w / 2, dec.within_state_component, w, color=ORANGE, label='Within-state')
    ax.scatter(xx, dec.observed_total, c='black', s=30, zorder=3, label='Observed total')
    ax.axhline(0, color='black', lw=0.8)
    ax.set(xticks=xx, xticklabels=dec.study_id, ylabel='Component contribution', title='Frozen stress-adjusted ECM coordinate, not the 25-gene score')
    ax.legend(frameon=False, ncol=3, loc='upper left')
    ax.grid(axis='y', alpha=0.16)
    ax.set_axisbelow(True)
    all_donors = read(AUDIT / 'All donor proportions.csv')
    groups = {'GSE251686': ('mild', 'severe'), 'GSE244889': ('MDD', 'SDD'), 'GSE153066': ('relatively_normal', 'degenerated')}
    for k, (study, (ref, dis)) in enumerate(groups.items()):
        ax = fig.add_subplot(grid[1, k])
        label(ax, 'BCD'[k])
        sub = all_donors[all_donors.study_id.eq(study)].copy()
        sub['disease'] = sub.condition.eq(dis)
        sub = sub.sort_values(['disease', 'remodeling_fraction'], ascending=[True, False]).reset_index(drop=True)
        yy = np.arange(len(sub))
        ax.scatter(sub.remodeling_fraction * 100, yy, c=[ORANGE if v else BLUE for v in sub.disease], s=28, zorder=3)
        for disease, col in [(False, BLUE), (True, ORANGE)]:
            ii = np.flatnonzero(sub.disease.eq(disease))
            mean = sub.loc[ii, 'remodeling_fraction'].mean() * 100
            ax.plot([mean, mean], [ii.min() - 0.35, ii.max() + 0.35], color=col, alpha=0.65, lw=1.4)
        ax.set(yticks=yy, yticklabels=sub.donor_id, title=study, xlim=(-2, 90), xlabel='State fraction (%)', ylim=(len(sub) - 0.4, -0.6))
        ax.grid(axis='x', alpha=0.16)
        ax.set_axisbelow(True)
    fig.legend(handles=[Line2D([], [], marker='o', linestyle='', color=BLUE, label='Reference / milder'), Line2D([], [], marker='o', linestyle='', color=ORANGE, label='Disease / more severe')], loc='lower center', bbox_to_anchor=(0.5, 0.03), ncol=2, frameon=False)
    fig.text(0.095, 0.012, 'B-D show every donor and group means from frozen fractions. The coordinate decomposition in A is a separate endpoint, not evidence that composition dominates.', fontsize=8)
    finish(fig, 'Supplementary Figure 22', SUPP)

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def writepdf(page, path):
    writer = PdfWriter()
    writer.add_page(page)
    with path.open('wb') as f:
        writer.write(f)

def pageoverlay(width, height, draw):
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(width, height), invariant=1)
    draw(c)
    c.save()
    return PdfReader(io.BytesIO(stream.getvalue())).pages[0]

def wrap_supplement(raw, number, destination):
    page = PdfReader(raw).pages[0]
    w, h = (float(page.mediabox.width), float(page.mediabox.height))
    if number <= 18:
        captions = json.loads((HERE / 'Inputs/Caption definitions.json').read_text(encoding='utf-8'))
        titles = json.loads((HERE / 'Inputs/Supplement titles.json').read_text(encoding='utf-8'))
        text = captions[f'S{number}']['caption']
        text = re.sub('^<b>.*?</b>\\s*', '', text, count=1) if text.startswith('<b>') and 'Supplementary Figure' in text[:160] else text
        para = Paragraph('<b>Figure legend.</b> ' + text, ParagraphStyle('caption', fontName='Helvetica', fontSize=11, leading=14.2, textColor='#111111'))
        _, ch = para.wrap(w, 2000)
        band = math.ceil(ch + 32)
        width, height = (w + 48, h + 72 + band)

        def draw(c):
            c.setFont('Helvetica-Bold', 11)
            c.drawString(24, height - 23, f'Supplementary Figure S{number}  {titles[number - 1]}')
            para.drawOn(c, 24, band + 12 - ch)
            c.setFont('Helvetica', 8)
            c.drawRightString(width - 24, 15, str(number))
        result = pageoverlay(width, height, draw)
        result.merge_transformed_page(page, Transformation().translate(tx=24, ty=32 + band))
    elif number <= 20:
        titles = {19: 'Literature marker context and cellular expression backgrounds', 20: 'Mixed phenotype prevalence and descriptive clinical context'}
        width, height = (w + 48, h + 72)

        def draw(c):
            c.setFont('Helvetica-Bold', 11)
            c.drawString(24, h + 48, f'Supplementary Figure S{number}  {titles[number]}')
            c.setFont('Helvetica', 8)
            c.drawString(24, 15, 'Previously completed expression-context audit.')
        result = pageoverlay(width, height, draw)
        result.merge_transformed_page(page, Transformation().translate(tx=24, ty=32))
    elif number == 21:
        width, height = (w + 48, h + 86)

        def draw(c):
            c.setFont('Helvetica-Bold', 11)
            c.drawString(24, h + 61, 'Supplementary Figure S21  Anatomical context of aligned ECM expression')
            c.setFont('Helvetica', 8)
            c.drawString(24, 16, 'Related anatomical expression does not establish equivalence of the complete fibrotic phenotype across tissues.')
        result = pageoverlay(width, height, draw)
        result.merge_transformed_page(page, Transformation().translate(tx=24, ty=32))
    else:
        titles = {22: 'ECM decomposition and donor state fractions', 23: 'Independent phenotype assessment in GSE165722', 24: 'Cluster 18 donor and detection depth sensitivity'}
        width, height = (w, h + 54)

        def draw(c):
            c.setFont('Helvetica-Bold', 11)
            c.drawString(24, height - 25, f'Supplementary Figure S{number}  {titles[number]}')
        result = pageoverlay(width, height, draw)
        result.merge_page(page)
    writepdf(result, destination)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Supplementary Figure S22 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S22.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        common_style()
        read = cached_read
        supplementary22()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 22, output)
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
