from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S15'
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
INPUT_FILES = {'cep-np/08_robustness/cell_level_coexpression_audit/donor_coexpression_proportions.tsv': 'Inputs/265d43dd8a donor coexpression proportions.tsv', 'cep-np/08_robustness/cell_level_coexpression_audit/dipper_target_protein_specimens.tsv': 'Inputs/b0004c5b82 dipper target protein specimens.tsv'}
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
import argparse, hashlib, json, re, shutil
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd
import seaborn as sns
ROOT = (HERE / 'Source keys/cep-np/scripts/61_polish_supplementary_figures.py').resolve().parents[1]
PKG = ROOT / 'manuscript_package_BMC_Genomics_v1'
OUT = PKG / 'figures/supplementary'
COEXP = ROOT / '08_robustness/cell_level_coexpression_audit'
ORDER = ['GSE251686', 'GSE153066', 'GSE244889', 'GSE233666', 'GSE222182', 'GSE160756_NP', 'GSE160756_AF', 'GSE160756_CEP']
EVENTS = ['double_positive', 'triple_A', 'triple_B']
EVENT_LABELS = {'double_positive': '$\\mathit{AEBP1}^{+}\\,\\mathit{COL5A2}^{+}$', 'triple_A': '$\\mathit{AEBP1}^{+}\\,\\mathit{COL5A2}^{+}\\,\\mathit{KRT19}^{+}$', 'triple_B': '$\\mathit{AEBP1}^{+}\\,\\mathit{COL5A2}^{+}\\,\\mathit{CA12}^{+}$'}
REPORT = []
SUFFIXES = {1: 'umap_GSE222182_primary', 2: 'umap_GSE222182_culture', 3: 'umap_GSE244889', 4: 'umap_GSE160756_tissue_background', 5: 'GSE160756_state_fraction_by_tissue', 6: 'GSE160756_frozen_signature_pseudobulk_heatmap', 7: 'umap_GSE233666_state_existence', 8: 'GSE186542_frozen_bulk_score', 9: 'state_alignment_heatmap', 10: 'independent_trajectory_scatter', 11: 'regulator_candidate_heatmap', 12: 'stratified_frozen_signature_meta_forest', 13: 'same_barcode_coexpression_effects', 14: 'ambient_proxy_and_UMI_threshold', 15: 'DIPPER_regional_protein_distribution', 16: 'leave_two_markers_out_audit', 17: 'detection_depth_standardization', 18: 'absolute_same_barcode_prevalence'}

def read(p):
    return pd.read_csv(p, sep='\t')

def cohort(s):
    return str(s).replace('GSE160756_', 'GSE160756 ').replace('_', ' ')

def style():
    sns.set_theme(style='ticks', context='paper')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'mathtext.fontset': 'dejavusans', 'savefig.dpi': 300})

def stem(n):
    return f'Supplementary_Figure_{n:02d}_{SUFFIXES[n]}'

def save(fig, n, note=''):
    name = stem(n)
    if note:
        fig.text(0.5, 0.012, note, ha='center', va='bottom', fontsize=8.5)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inactive_tick_text = set()
    for ax in fig.axes:
        for axis, limits in [(ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())]:
            low, high = sorted(limits)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not low - 1e-08 <= tick.get_loc() <= high + 1e-08:
                    inactive_tick_text.update([id(tick.label1), id(tick.label2)])
    outside = []
    for t in fig.findobj(Text):
        if id(t) in inactive_tick_text or not t.get_visible() or (not t.get_text().strip()):
            continue
        b = t.get_window_extent(renderer)
        if b.x0 < -2 or b.y0 < -2 or b.x1 > fig.bbox.width + 2 or (b.y1 > fig.bbox.height + 2):
            outside.append(t.get_text())
    fig.savefig(OUT / f'{name}.png', dpi=300, facecolor='white')
    fig.savefig(OUT / f'{name}.pdf', facecolor='white')
    REPORT.append({'figure': n, 'stem': name, 'canvas_inches': list(fig.get_size_inches()), 'out_of_canvas_text': outside, 'presentation_only': True})
    plt.close(fig)

def ambient_and_dipper():
    df = read(COEXP / 'donor_coexpression_proportions.tsv')
    df['analysis_cohort'] = np.where(df.study.eq('GSE160756'), df.study + '_' + df.tissue, df.study)
    pops = ['outside_broad_mask', 'nonremodeling_NP_related', 'remodeling']
    agg = df[df.event.isin(EVENTS) & df.population.isin(pops)].groupby(['analysis_cohort', 'threshold_umi', 'event', 'population'], observed=True).proportion.median().reset_index()
    colors = dict(zip(ORDER, sns.color_palette('colorblind', len(ORDER))))
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 8))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.19, top=0.85, hspace=0.43, wspace=0.3)
    for r, threshold in enumerate([1, 2]):
        for c, event in enumerate(EVENTS):
            ax = axes[r, c]
            block = agg[agg.event.eq(event) & agg.threshold_umi.eq(threshold)]
            for k, co in enumerate(ORDER):
                b = block[block.analysis_cohort.eq(co)].set_index('population').reindex(pops)
                ax.scatter(np.arange(3) + (k - 3.5) * 0.045, b.proportion, color=colors[co], s=23)
            ax.set_xticks(range(3), labels=['Outside broad\ncompartment', 'Other NP-\nrelated cells', 'Remodelling\nstate'])
            ax.set_title(EVENT_LABELS[event] + '\n' + '$\\geq$' + f' {threshold} UMI per gene', fontsize=10)
            ax.set_ylabel('Median within-stratum fraction' if c == 0 else '')
    handles = [Line2D([], [], marker='o', linestyle='', color=colors[k], label=cohort(k), markersize=4) for k in ORDER]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.045), ncol=4, frameon=False, fontsize=8.5)
    fig.suptitle('Ambient-RNA proxy and UMI-threshold sensitivity', fontsize=13, fontweight='bold', y=0.97)
    save(fig, 14)
    dp = read(COEXP / 'dipper_target_protein_specimens.tsv')
    fig, axes = plt.subplots(1, 3, figsize=(11.3, 5.2))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.8, bottom=0.23, wspace=0.35)
    for i, (ax, gene) in enumerate(zip(axes, ['AEBP1', 'COL5A2', 'KRT19'])):
        b = dp[dp.gene.eq(gene)]
        regions = ['OAF', 'IAF', 'IAF_NP', 'NP']
        sns.boxplot(data=b, x='compartment', y='log2_abundance', hue='age_group', order=regions, hue_order=['young', 'old'], palette=['#4C78A8', '#F2A269'], whis=(0, 100), fliersize=0, ax=ax)
        for j, region in enumerate(regions):
            for k, age in enumerate(['young', 'old']):
                block = b[b.compartment.eq(region) & b.age_group.eq(age)]
                ax.scatter(j + (k - 0.5) * 0.38 + np.linspace(-0.075, 0.075, len(block)), block.log2_abundance, s=9, color='#666666', alpha=0.55, zorder=3)
        ax.set_xticks(range(4), labels=['OAF', 'IAF', 'NP/IAF', 'NP'])
        ax.set(xlabel='DIPPER region', ylabel='Log2 protein abundance', title=gene + ('*' if gene == 'KRT19' else ''))
        if i == 0:
            ax.legend(title='Age group', frameon=False, fontsize=8)
        else:
            ax.get_legend().remove()
    fig.suptitle('Regional protein distribution in DIPPER', fontsize=13, fontweight='bold', y=0.95)
    save(fig, 15, '* KRT19 source protein group is flagged as a contaminant. Regional abundance does not establish cellular colocalization.')

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
    raw = WORK_DIR / 'Supplementary Figure S15 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S15.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        style()

        def selected_save(fig, n, note=''):
            if n == 15:
                if note:
                    fig.text(0.5, 0.012, note, ha='center', va='bottom', fontsize=8.5)
                rawsave(fig, raw, facecolor='white')
            plt.close(fig)
        save = selected_save
        ambient_and_dipper()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 15, output)
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
