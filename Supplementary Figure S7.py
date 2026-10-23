from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S7'
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
INPUT_FILES = {'cep-np/10_figures_manuscript/supplementary_polish_inputs/s07_frozen_umap_display.tsv.gz': 'Inputs/2844493764 S07 frozen umap display.tsv.gz', 'cep-np/10_figures_manuscript/supplementary_polish_inputs/display_input_provenance.json': 'Inputs/e450057c2e display input provenance.json'}
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
INPUT = ROOT / '10_figures_manuscript/supplementary_polish_inputs'
REPORT = []
SUFFIXES = {1: 'umap_GSE222182_primary', 2: 'umap_GSE222182_culture', 3: 'umap_GSE244889', 4: 'umap_GSE160756_tissue_background', 5: 'GSE160756_state_fraction_by_tissue', 6: 'GSE160756_frozen_signature_pseudobulk_heatmap', 7: 'umap_GSE233666_state_existence', 8: 'GSE186542_frozen_bulk_score', 9: 'state_alignment_heatmap', 10: 'independent_trajectory_scatter', 11: 'regulator_candidate_heatmap', 12: 'stratified_frozen_signature_meta_forest', 13: 'same_barcode_coexpression_effects', 14: 'ambient_proxy_and_UMI_threshold', 15: 'DIPPER_regional_protein_distribution', 16: 'leave_two_markers_out_audit', 17: 'detection_depth_standardization', 18: 'absolute_same_barcode_prevalence'}
UMAPS = [(1, 'GSE222182/GSE222182_primary_independent_annotated.h5ad', 'GSE222182 primary tissue', [('sample_id', 'Sample'), ('condition', 'Symptom context'), ('leiden_0.8', 'Independent cluster'), ('transition_score', 'Frozen remodelling score'), ('broad_NP_mesenchymal', 'Broad NP/mesenchymal compartment'), ('eligible_remodeling_state', 'Remodelling state')]), (2, 'GSE222182/GSE222182_culture_independent_annotated.h5ad', 'GSE222182 cultured cells', [('sample_id', 'Sample'), ('condition', 'Culture condition'), ('leiden_0.8', 'Independent cluster'), ('transition_score', 'Frozen remodelling score'), ('broad_NP_mesenchymal', 'Broad NP/mesenchymal compartment'), ('eligible_remodeling_state', 'Remodelling state')]), (3, 'GSE244889/GSE244889_independent_annotated.h5ad', 'GSE244889', [('donor_id', 'Donor'), ('condition', 'Degeneration group'), ('leiden_0.8', 'Independent cluster'), ('transition_score', 'Frozen remodelling score'), ('broad_NP_mesenchymal', 'Broad NP/mesenchymal compartment'), ('eligible_remodeling_state', 'Remodelling state')]), (4, 'GSE160756/GSE160756_tissue_background_annotated.h5ad', 'GSE160756', [('sample_id', 'Library'), ('tissue', 'Tissue'), ('leiden_unintegrated', 'Independent cluster (unintegrated)'), ('transition_score', 'Frozen remodelling score'), ('remodeling_state_unintegrated', 'Remodelling state (unintegrated)')]), (7, 'GSE233666/GSE233666_independent_state_existence.h5ad', 'GSE233666', [('sample_id', 'Patient'), ('age', 'Age (years)'), ('leiden_0.8', 'Independent cluster'), ('transition_score', 'Frozen remodelling score'), ('remodeling_state_main', 'Remodelling state')])]

def label(s):
    s = str(s).replace('_', ' ')
    if s.startswith(('asymptomatic ', 'back pain ', 'no stress ', 'stress ')):
        s = s[0].upper() + s[1:]
    return {'MDD': 'Mild degeneration', 'SDD': 'Severe degeneration', 'asymptomatic': 'Asymptomatic', 'back pain': 'Back pain', 'no stress': 'No stress', 'stress': 'Stress', 'early': 'Early', 'advanced': 'Advanced'}.get(s, s)

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

def umaps():
    for n, path, title, fields in UMAPS:
        df = pd.read_csv(INPUT / f'S{n:02d}_frozen_umap_display.tsv.gz', sep='\t', float_precision='round_trip')
        fig = plt.figure(figsize=(10.8, 13.8))
        gs = fig.add_gridspec(3, 2, left=0.075, right=0.94, bottom=0.035, top=0.95, hspace=0.21, wspace=0.24)
        for i, (field, heading) in enumerate(fields):
            inner = gs[i // 2, i % 2].subgridspec(2, 1, height_ratios=[3.1, 0.82], hspace=0.16)
            ax = fig.add_subplot(inner[0])
            la = fig.add_subplot(inner[1])
            la.axis('off')
            ax.set_title(heading, pad=9)
            ax.text(-0.12, 1.06, chr(65 + i), transform=ax.transAxes, fontweight='bold', fontsize=13)
            ax.set_xlabel('UMAP 1', labelpad=3)
            ax.set_ylabel('UMAP 2', labelpad=3)
            ax.set_xticks([])
            ax.set_yticks([])
            vals = df[field]
            if field in ['age', 'transition_score']:
                drawing_order = np.argsort(vals.to_numpy(), kind='stable')
                points = ax.scatter(df.UMAP1.iloc[drawing_order], df.UMAP2.iloc[drawing_order], c=vals.iloc[drawing_order], s=1.8, linewidths=0, cmap='viridis', rasterized=True)
                cb = fig.colorbar(points, ax=ax, fraction=0.035, pad=0.02)
                cb.ax.tick_params(labelsize=8)
                cb.outline.set_linewidth(0.5)
            else:
                levels = list(dict.fromkeys(vals.astype(str)))
                if field.startswith('leiden'):
                    levels = sorted(levels, key=lambda v: int(v))
                elif vals.dtype == bool:
                    levels = ['False', 'True'] if vals.any() and (~vals).any() else sorted(levels)
                else:
                    levels = sorted(levels)
                palette = sns.color_palette('tab20', len(levels)) if len(levels) > 10 else sns.color_palette('tab10', len(levels))
                if len(levels) > 20:
                    palette = sns.color_palette('husl', len(levels))
                if vals.dtype == bool:
                    palette = [{'False': '#3182BD', 'True': '#F28E2B'}[v] for v in levels]
                color_map = dict(zip(levels, palette))
                point_colors = [color_map[v] for v in vals.astype(str)]
                ax.scatter(df.UMAP1, df.UMAP2, c=point_colors, s=1.8, linewidths=0, rasterized=True)
                handles = []
                for k, color in zip(levels, palette):
                    text = {'False': 'No', 'True': 'Yes'}.get(k, label(k))
                    handles.append(Line2D([], [], marker='o', linestyle='', markersize=4, color=color, label=text))
                cols = 7 if field.startswith('leiden') else 2 if max(map(lambda h: len(h.get_label()), handles)) > 10 else min(4, len(handles))
                la.legend(handles=handles, loc='upper center', ncol=cols, frameon=False, fontsize=8.3, columnspacing=0.9, handletextpad=0.4, borderaxespad=0)
            sns.despine(ax=ax)
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.985)
        save(fig, n)

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
    raw = WORK_DIR / 'Supplementary Figure S7 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S7.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        style()

        def selected_save(fig, n, note=''):
            if n == 7:
                if note:
                    fig.text(0.5, 0.012, note, ha='center', va='bottom', fontsize=8.5)
                rawsave(fig, raw, facecolor='white')
            plt.close(fig)
        save = selected_save
        UMAPS = [x for x in UMAPS if x[0] == 7]
        umaps()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 7, output)
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
