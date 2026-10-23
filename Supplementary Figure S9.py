from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S9'
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
INPUT_FILES = {'cep-np/08_robustness/state_deep_dive/outcome_blind_cluster_state_alignment.tsv': 'Inputs/b55f668b08 outcome blind cluster state alignment.tsv', 'cep-np/08_robustness/state_deep_dive/independent_diffusion_pseudotime_clusters.tsv': 'Inputs/e1e3ef66a3 independent diffusion pseudotime clusters.tsv', 'cep-np/08_robustness/state_deep_dive/exploratory_regulator_priority.tsv': 'Inputs/2b4dc31532 exploratory regulator priority.tsv'}
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
DEEP = ROOT / '08_robustness/state_deep_dive'
REPORT = []
SUFFIXES = {1: 'umap_GSE222182_primary', 2: 'umap_GSE222182_culture', 3: 'umap_GSE244889', 4: 'umap_GSE160756_tissue_background', 5: 'GSE160756_state_fraction_by_tissue', 6: 'GSE160756_frozen_signature_pseudobulk_heatmap', 7: 'umap_GSE233666_state_existence', 8: 'GSE186542_frozen_bulk_score', 9: 'state_alignment_heatmap', 10: 'independent_trajectory_scatter', 11: 'regulator_candidate_heatmap', 12: 'stratified_frozen_signature_meta_forest', 13: 'same_barcode_coexpression_effects', 14: 'ambient_proxy_and_UMI_threshold', 15: 'DIPPER_regional_protein_distribution', 16: 'leave_two_markers_out_audit', 17: 'detection_depth_standardization', 18: 'absolute_same_barcode_prevalence'}

def read(p):
    return pd.read_csv(p, sep='\t')

def gene_text(s):
    out = str(s).replace('remodeling', 'remodelling')
    pattern = '(?<![A-Za-z0-9])([A-Z][A-Z0-9]{1,})(?![A-Za-z0-9])'

    def convert(m):
        token = m[1]
        if token in ['NP', 'AF', 'CEP', 'ECM', 'UMI', 'RNA', 'DNA', 'DIPPER', 'COL1', 'COL3'] or token.startswith('GSE') or re.fullmatch('C\\d+', token):
            return token
        return '$\\mathit{' + token + '}$'
    return re.sub(pattern, convert, out)

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

def deep_dive():
    df = read(DEEP / 'outcome_blind_cluster_state_alignment.tsv')
    cols = [c for c in df if c.startswith('score_')]
    names = [gene_text(f'{r.study_id} C{r.cluster} | {r.state_detail_original}') for r in df.itertuples()]
    xlabels = [gene_text(c.removeprefix('score_').replace('_', ' ').replace('np like', 'NP-like').replace('mgp', 'MGP').replace('msmp', 'MSMP').replace('chi3l1', 'CHI3L1').replace('chi3l2', 'CHI3L2').replace('prg4', 'PRG4').replace('cilp', 'CILP').replace('ecm', 'ECM')) for c in cols]
    fig, ax = plt.subplots(figsize=(12.8, 14.5))
    fig.subplots_adjust(left=0.53, right=0.97, bottom=0.19, top=0.95)
    sns.heatmap(df[cols], cmap='vlag', center=0, ax=ax, yticklabels=names, xticklabels=xlabels, cbar_kws={'label': 'Within-cohort signature score', 'shrink': 0.65})
    ax.tick_params(axis='y', labelsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=55, ha='right', fontsize=8)
    ax.set(xlabel='Frozen discovery state', ylabel='', title='Positive-marker alignment across NP-related clusters')
    save(fig, 9)
    tr = read(DEEP / 'independent_diffusion_pseudotime_clusters.tsv')
    studies = ['GSE251686', 'GSE229711', 'GSE230808', 'GSE153066']
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 10.2))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.9, bottom=0.14, hspace=0.38, wspace=0.3)
    from adjustText import adjust_text
    for i, (ax, study) in enumerate(zip(axes.flat, studies)):
        b = tr[tr.study_id.eq(study)].copy()
        sizes = 40 + (b.cells - b.cells.min()) / max(b.cells.max() - b.cells.min(), 1) * 240
        target = b.eligible_remodeling.astype(str).str.lower().eq('true')
        ax.scatter(b.median_pseudotime, b.target_signature_score, s=sizes, c=np.where(target, '#D7301F', '#3182BD'), edgecolors='white', linewidths=0.5)
        ax.margins(x=0.17, y=0.2)
        if study == 'GSE230808':
            positions = {1: (0.01, 0.43), 12: (0.15, 0.21), 4: (-0.055, 0.15), 2: (0.15, 0.11), 0: (0.16, 0.035), 11: (-0.055, -0.02), 10: (0.17, -0.045), 20: (0.12, -0.11), 5: (-0.04, -0.25)}
            for _, row in b.iterrows():
                cid = int(row['leiden_0.8'])
                xy = (row.median_pseudotime, row.target_signature_score)
                target_pos = positions.get(cid, (xy[0] + 0.015, xy[1] + 0.025))
                ax.annotate(str(cid), xy, xytext=target_pos, fontsize=9, ha='center', va='center', arrowprops=dict(arrowstyle='-', color='#777777', lw=0.55, shrinkA=2, shrinkB=6))
        else:
            texts = [ax.text(row.median_pseudotime, row.target_signature_score, str(int(row['leiden_0.8'])), fontsize=9) for _, row in b.iterrows()]
            adjust_text(texts, ax=ax, expand=(1.35, 1.6), force_text=(0.7, 0.9), force_static=(0.5, 0.7), iter_lim=400, ensure_inside_axes=True, arrowprops=dict(arrowstyle='-', color='#777777', lw=0.5))
        ax.set(xlabel='Median diffusion pseudotime', ylabel='Remodelling signature score', title=study)
        ax.text(-0.13, 1.06, chr(65 + i), transform=ax.transAxes, fontweight='bold', fontsize=13)
    handles = [Line2D([], [], marker='o', linestyle='', color=c, label=t) for c, t in [('#D7301F', 'Eligible remodelling cluster'), ('#3182BD', 'Other NP-related cluster')]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.055), ncol=2, frameon=False)
    fig.suptitle('Independent diffusion-pseudotime comparisons', fontsize=14, fontweight='bold', y=0.975)
    save(fig, 10, 'Numbers are frozen cluster IDs; point size reflects cell count. Pseudotime is not a disease timeline.')
    reg = read(DEEP / 'exploratory_regulator_priority.tsv').head(12)
    values = reg.set_index('gene_symbol')[['discovery_target_minus_other', 'GSE230808_target_minus_other', 'GSE153066_target_minus_other']]
    values.columns = ['GSE251686', 'GSE230808', 'GSE153066']
    fig, ax = plt.subplots(figsize=(7.4, 8.2))
    fig.subplots_adjust(left=0.2, right=0.94, bottom=0.13, top=0.91)
    annotations = values.map(lambda v: '0.00' if abs(v) < 0.005 else f'{v:.2f}')
    sns.heatmap(values, ax=ax, cmap='vlag', center=0, annot=annotations, fmt='', cbar_kws={'label': 'Target minus other log-expression', 'shrink': 0.8})
    ax.set_yticklabels(values.index, rotation=0, fontstyle='italic')
    ax.set_xticklabels(values.columns, rotation=0)
    ax.set(xlabel='Cohort', ylabel='Candidate regulator gene', title='Exploratory regulator expression contrasts')
    save(fig, 11)

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
    raw = WORK_DIR / 'Supplementary Figure S9 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S9.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        style()

        def selected_save(fig, n, note=''):
            if n == 9:
                if note:
                    fig.text(0.5, 0.012, note, ha='center', va='bottom', fontsize=8.5)
                rawsave(fig, raw, facecolor='white')
            plt.close(fig)
        save = selected_save
        deep_dive()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 9, output)
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
