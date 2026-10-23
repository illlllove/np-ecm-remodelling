from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S21'
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
INPUT_FILES = {'cep-np/manuscript_package_bmc_genomics_v3/editorial/presentation inputs/cluster markers.tsv.gz': 'Inputs/0613a079f0 cluster markers.tsv.gz', 'cep-np/07_validation/gse233666/cluster_markers.tsv.gz': 'Inputs/0613a079f0 cluster markers.tsv.gz', 'cep-np/07_validation/gse160756/unintegrated_sample_tissue_abundance.tsv': 'Inputs/e889a13764 unintegrated sample tissue abundance.tsv', 'cep-np/07_validation/gse160756/unintegrated_tissue_summary.tsv': 'Inputs/4519039a01 unintegrated tissue summary.tsv', 'cep-np/07_validation/gse160756/harmony_tissue_summary.tsv': 'Inputs/5bd0f7afa7 harmony tissue summary.tsv', 'cep-np/07_validation/gse233666/sample_state_presence.tsv': 'Inputs/49daf54490 sample state presence.tsv'}
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

def figure2():
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.8), gridspec_kw={'height_ratios': [1.05, 0.95]})
    ax = axes[0, 0]
    label_panel(ax, 'A')
    markers = pd.read_csv(ROOT / '07_validation/GSE233666/cluster_markers.tsv.gz', sep='\t')
    genes = ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2', 'COL6A1', 'COL6A2', 'SPARC', 'POSTN', 'PCOLCE', 'MMP2', 'MMP14', 'CTHRC1', 'FAP', 'AEBP1']
    q = markers.loc[markers.group.astype(str).isin(['5', '11']) & markers.names.isin(genes)].copy()
    pivot = q.pivot_table(index='names', columns='group', values='logfoldchanges', aggfunc='first').reindex(genes).fillna(0)
    pivot.columns = [f'State {c}' for c in pivot.columns]
    sns.heatmap(pivot, cmap='YlOrRd', linewidths=0.4, linecolor='white', cbar_kws={'label': 'cluster log fold-change'}, ax=ax)
    ax.set_title('GSE233666 state-defining matrix programme')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax = axes[0, 1]
    label_panel(ax, 'B')
    presence = pd.read_csv(ROOT / '07_validation/GSE233666/sample_state_presence.tsv', sep='\t')
    presence = presence.sort_values('main_state_fraction')
    ax.barh(presence.sample_id, presence.main_state_fraction, color=BLUE, alpha=0.85)
    for i, row in enumerate(presence.itertuples()):
        ax.text(row.main_state_fraction + 0.006, i, f'{row.main_state_fraction:.3f}', va='center', fontsize=8)
    ax.set_xlim(0, max(presence.main_state_fraction) * 1.25)
    ax.set_xlabel('state fraction within broad NP/mesenchymal cells')
    ax.set_title('State detected in all four herniated-NP patients')
    ax = axes[1, 0]
    label_panel(ax, 'C')
    tissue = pd.read_csv(ROOT / '07_validation/GSE160756/unintegrated_sample_tissue_abundance.tsv', sep='\t')
    order = ['NP', 'AF', 'CEP']
    sns.stripplot(data=tissue, x='tissue', y='remodeling_fraction', order=order, hue='sample_id', size=8, jitter=0.08, ax=ax)
    med = tissue.groupby('tissue').remodeling_fraction.median().reindex(order)
    for i, v in enumerate(med):
        ax.plot([i - 0.22, i + 0.22], [v, v], color='black', lw=2)
    ax.get_legend().remove()
    ax.set_xlabel('')
    ax.set_ylabel('primary state fraction')
    ax.set_title('GSE160756 primary tissue-background analysis')
    ax = axes[1, 1]
    label_panel(ax, 'D')
    rows = []
    for rep in ['unintegrated', 'harmony']:
        s = pd.read_csv(ROOT / f'07_validation/GSE160756/{rep}_tissue_summary.tsv', sep='\t')
        for r in s.itertuples():
            rows.append({'representation': 'Primary' if rep == 'unintegrated' else 'Harmony sensitivity', 'tissue': r.tissue, 'median': r.median_fraction})
    d = pd.DataFrame(rows)
    sns.barplot(data=d, x='tissue', y='median', hue='representation', order=order, palette=[BLUE, ORANGE], ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel('median state fraction')
    ax.set_title('NP enrichment persists across representations')
    ax.legend(frameon=False, title='')
    fig.suptitle('Independent data identify an NP-enriched activated fibrocartilage-like ECM-remodelling state', fontsize=14, fontweight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, 'Figure2_independent_identity_and_tissue_background')

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
    raw = WORK_DIR / 'Supplementary Figure S21 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S21.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        style()

        def local_save(fig, name, output_dir=None):
            rawsave(fig, raw, bbox_inches='tight', facecolor='white')
            plt.close(fig)
        save = local_save
        figure2()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    cropped = WORK_DIR / 'S21 journal crop.pdf'
    crop(raw, cropped)
    wrap_supplement(cropped, 21, output)
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
