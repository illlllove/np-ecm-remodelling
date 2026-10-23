from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S24'
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
INPUT_FILES = {'results/gse165722_independent_validation/06_cluster18_sensitivity/tables/paired_donor_effects.tsv': 'Inputs/576a5f3d7c paired donor effects.tsv', 'results/gse165722_independent_validation/06_cluster18_sensitivity/tables/depth_standardized_by_donor.tsv': 'Inputs/d634683d89 depth standardized by donor.tsv', 'results/gse165722_independent_validation/06_cluster18_sensitivity/tables/donor_inference_summary.tsv': 'Inputs/e8b7ca0360 donor inference summary.tsv', 'results/gse165722_independent_validation/06_cluster18_sensitivity/tables/leave_one_donor_out_effects.tsv': 'Inputs/0978a9e1fe leave one donor out effects.tsv', 'results/gse165722_independent_validation/06_cluster18_sensitivity/tables/absolute_prevalence.tsv': 'Inputs/6258ab4b24 absolute prevalence.tsv'}
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
G165 = HOMEBASE / 'results/GSE165722_independent_validation'
MAIN = PKG / 'figures/main'
SUPP = PKG / 'figures/supplementary'
SOURCES = {}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def picture(path):
    path = Path(path)
    SOURCES[str(path)] = sha(path)
    return plt.imread(path)

def finish(fig, name, folder=MAIN):
    for ext in ['pdf', 'png']:
        fig.savefig(folder / f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def common_style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10.5, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

def supplementary24():
    source = G165 / '06_cluster18_sensitivity/figures/cluster18_donor_and_depth_audit.png'
    fig, ax = plt.subplots(figsize=(10, 9))
    fig.subplots_adjust(left=0.025, right=0.99, top=0.94, bottom=0.12)
    ax.imshow(picture(source))
    ax.set_axis_off()
    fig.suptitle('Cluster 18 focused post-selection descriptive sensitivity', fontsize=12, y=0.985)
    fig.text(0.035, 0.074, 'The original cluster 18 graphic is unchanged. Paired co-detection: five donors, nominal P=0.0625.', fontsize=8)
    fig.text(0.035, 0.048, 'Depth common support: two donors; retained effects 42.4% and 34.9%, nominal P=0.5. Depth sensitivity is material.', fontsize=8)
    fig.text(0.035, 0.023, 'Full frozen identity LODO passed in 7/8 omissions. This is not the union endpoint, a clean validation pass or a disease-effect test.', fontsize=8)
    finish(fig, 'Supplementary Figure 24', SUPP)

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

def cluster18_plot_block():
    from pathlib import Path
    import pandas as pd, numpy as np, json, shutil, re, datetime
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ROOT = HERE / 'Source keys/results/GSE165722_independent_validation'
    O = ROOT / '06_cluster18_sensitivity'
    T = O / 'tables'

    def rd(n):
        return pd.read_csv(T / n, sep='\t')
    p = rd('paired_donor_effects.tsv')
    de = rd('depth_standardized_by_donor.tsv')
    s = rd('donor_inference_summary.tsv')
    lo = rd('leave_one_donor_out_effects.tsv')
    ab = rd('absolute_prevalence.tsv')
    q = p[(p.contrast == 'original_non_target_NP') & (p.count_threshold == 1)]
    eligible = q[q.paired_eligible]
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.4), gridspec_kw={'height_ratios': [1.5, 1]})
    ax = axes[0]
    x = np.arange(len(eligible))
    labels = []
    for j, (_, r) in enumerate(eligible.iterrows()):
        ax.plot([j, j], [r.comparator_fraction * 100, r.cluster18_fraction * 100], color='#94A3B8', lw=1)
        labels.append(r.donor[-3:] + '\n' + str(r.n_cluster18) + ' / ' + str(r.n_comparator))
    ax.scatter(x, eligible.cluster18_fraction * 100, s=50, color='#0F766E', label='Cluster 18')
    ax.scatter(x, eligible.comparator_fraction * 100, s=40, color='#475569', marker='s', label='Frozen non-target NP')
    ax.set_xticks(x, labels, fontsize=9)
    ax.set_ylabel('Double-positive cells (%)', fontsize=10)
    ax.set_xlabel('GSM suffix; cluster 18 / comparison cells', fontsize=9)
    ax.set_title('A  Within-donor co-detection (supplied counts >=1 per gene)', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc='upper left')
    ax.set_ylim(0, 100)
    ax = axes[1]
    v = s[(s.contrast == 'original_non_target_NP') & (s.analysis == 'depth_standardized')].sort_values('count_threshold')
    ax.bar([0, 1], v.effect_retained_ratio * 100, color=['#B45309', '#92400E'], width=0.55)
    ax.axhline(50, color='#991B1B', ls='--', lw=1, label='Frozen 50% retention boundary')
    for i, r in enumerate(v.itertuples()):
        ax.text(i, r.effect_retained_ratio * 100 + 3, f'{r.effect_retained_ratio:.1%}', ha='center', fontsize=10)
    ax.set_xticks([0, 1], ['>=1 count each', '>=2 counts each'])
    ax.set_ylabel('Effect retained (%)')
    ax.set_title('B  Depth-standardized / matched unadjusted donor-mean RD (n=2)', loc='left', fontsize=11)
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(h_pad=2.5)
    fig.savefig(O / 'figures/cluster18_donor_and_depth_audit.png', dpi=220)
    fig.savefig(O / 'figures/cluster18_donor_and_depth_audit.svg')
    plt.close(fig)

def nested_cluster18():
    output = WORK_DIR / 'Cluster 18 rebuilt audit.png'
    old_save = Figure.savefig

    def save_local(fig, path, *args, **kwargs):
        if str(path).endswith('.png'):
            rawsave(fig, output, **kwargs)
    Figure.savefig = save_local
    try:
        cluster18_plot_block()
    finally:
        Figure.savefig = old_save
    return output

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Supplementary Figure S24 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S24.pdf'
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
        nested = nested_cluster18()
        plt.rcParams.update(prior)
        Figure.savefig = onlysave
        picture = lambda ignored: plt.imread(nested)
        supplementary24()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 24, output)
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
