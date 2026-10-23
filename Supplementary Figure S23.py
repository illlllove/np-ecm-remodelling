from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S23'
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
INPUT_FILES = {'results/gse165722_independent_validation/04_phenotype/same_barcode_events_by_donor_population.tsv': 'Inputs/fd4c70d908 same barcode events by donor population.tsv'}
MANIFEST = {'canonicalRoot': str(HERE / 'Source keys/cep-np/manuscript package BMC Genomics v6'), 'GSE165722display': {'cache': 'Inputs/GSE165722 frozen UMAP display.tsv.gz', 'sha256': '59efa25f08e8a7a6091b683c60d01a49b957a28cbf8de6c67b950a072e4f96d9', 'cells': 41694}}

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
G165 = HOMEBASE / 'results/GSE165722_independent_validation'
MAIN = PKG / 'figures/main'
SUPP = PKG / 'figures/supplementary'
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

def picture(path):
    path = Path(path)
    SOURCES[str(path)] = sha(path)
    return plt.imread(path)

def label(ax, letter, x=-0.13, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=15, fontweight='bold', ha='left', va='bottom')

def finish(fig, name, folder=MAIN):
    for ext in ['pdf', 'png']:
        fig.savefig(folder / f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def common_style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10.5, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

def supplementary23():
    fig = plt.figure(figsize=(14.8, 10))
    grid = fig.add_gridspec(2, 1, left=0.08, right=0.97, bottom=0.12, top=0.94, height_ratios=[1, 1.35], hspace=0.28)
    a = fig.add_subplot(grid[0])
    label(a, 'A', x=-0.04)
    a.imshow(picture(G165 / '05_figures/independent_umap_clusters_identity_donor.png'))
    a.set_axis_off()
    a.set_title('Independent clustering and cellular backgrounds in GSE165722', pad=12)
    b = fig.add_subplot(grid[1])
    label(b, 'B', x=-0.06)
    events = read(G165 / '04_phenotype/same_barcode_events_by_donor_population.tsv', sep='\t')
    q = events[events.event.eq('double_positive') & events.count_threshold.eq(1)]
    donors = sorted(q.donor.unique())
    colors = {'target': BLUE, 'non_target_NP': ORANGE, 'other_cells': GREEN}
    names = {'target': 'Frozen union: clusters 15 / 16 / 18', 'non_target_NP': 'Frozen non-target NP', 'other_cells': 'Outside NP mask'}
    for k, pop in enumerate(['target', 'non_target_NP', 'other_cells']):
        z = q[q.population.eq(pop)].set_index('donor').reindex(donors)
        b.scatter(np.arange(len(donors)) + (k - 1) * 0.13, z.fraction * 100, c=colors[pop], s=35, label=names[pop], zorder=3)
    b.set(xticks=range(len(donors)), xticklabels=[s.replace('GSM5048', '') for s in donors], ylabel='Double-positive cells (%)', xlabel='GSM suffix; donor/library identity only', ylim=(-2, 100), title='$\\mathit{AEBP1}$ + $\\mathit{COL5A2}$ in supplied counts at least 1 per gene')
    b.legend(frameon=False, ncol=3, loc='upper left', fontsize=8)
    b.grid(axis='y', alpha=0.16)
    b.set_axisbelow(True)
    fig.text(0.08, 0.068, 'Union endpoint only: six evaluable paired donors, nominal P=0.03125; depth-standardized five donors, P=0.0625.', fontsize=9)
    fig.text(0.08, 0.041, 'These P values do not apply to cluster 18 or a purified NP subset. Candidate identities are heterogeneous; clinical grade mapping is unresolved.', fontsize=8)
    fig.text(0.08, 0.018, 'Conditional supplied-count phenotype audit: no grade association, full 25-gene score or stable-discrete-subtype validation is asserted.', fontsize=8)
    finish(fig, 'Supplementary Figure 23', SUPP)

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

def nested_umap():
    table = ORIGINALREAD(HERE / MANIFEST['GSE165722display']['cache'], sep='\t', float_precision='round_trip', dtype={'leiden_0.8': str})
    palettes = json.loads((HERE / 'Inputs/GSE165722 palettes.json').read_text(encoding='utf-8'))['palettes']
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'savefig.bbox': 'tight'})
    import anndata as ad
    import scanpy as sc
    columns = ['leiden_0.8', 'independent_annotation_leiden_0.8', 'audit_donor']
    obj = ad.AnnData(X=np.empty((len(table), 0)), obs=table[columns].copy())
    obj.obsm['X_umap'] = table[['UMAP1', 'UMAP2']].to_numpy()
    for column, palette in zip(columns, palettes):
        obj.obs[column] = pd.Categorical(obj.obs[column].astype(str), categories=list(palette))
        obj.uns[column + '_colors'] = list(palette.values())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, column, title in zip(axes, columns, ['Independent clusters', 'Frozen phenotype annotation', 'Donor']):
        sc.pl.umap(obj, color=column, ax=ax, show=False, title=title, size=4, frameon=False, legend_loc='right margin')
    fig.tight_layout()
    target = WORK_DIR / 'GSE165722 rebuilt UMAP.png'
    rawsave(fig, target, dpi=200)
    plt.close(fig)
    return target

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Supplementary Figure S23 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S23.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        common_style()
        read = cached_read
        Figure.savefig = old_save
        prior = plt.rcParams.copy()
        plt.rcdefaults()
        nested = nested_umap()
        plt.rcParams.update(prior)
        Figure.savefig = onlysave
        picture = lambda ignored: plt.imread(nested)
        supplementary23()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 23, output)
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
