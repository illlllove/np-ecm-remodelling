from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S18'
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
INPUT_FILES = {'cep-np/08_robustness/prepublication_coexpression_qc/depth_standardized_cohort_summary.tsv': 'Inputs/3bff1a9435 depth standardized cohort summary.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/absolute_prevalence_cohort_summary.tsv': 'Inputs/0b6331f41f absolute prevalence cohort summary.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/leave_two_out_cohort_summary.tsv': 'Inputs/1bbb6f02d1 leave two out cohort summary.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/leave_two_out_state_support.tsv': 'Inputs/dc40c112a5 leave two out state support.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/depth_metrics_cohort_summary.tsv': 'Inputs/083963df8a depth metrics cohort summary.tsv'}
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
QC = ROOT / '08_robustness/prepublication_coexpression_qc'
ORDER = ['GSE251686', 'GSE153066', 'GSE244889', 'GSE233666', 'GSE222182', 'GSE160756_NP', 'GSE160756_AF', 'GSE160756_CEP']
EVENTS = ['double_positive', 'triple_A', 'triple_B']
EVENT_LABELS = {'double_positive': '$\\mathit{AEBP1}^{+}\\,\\mathit{COL5A2}^{+}$', 'triple_A': '$\\mathit{AEBP1}^{+}\\,\\mathit{COL5A2}^{+}\\,\\mathit{KRT19}^{+}$', 'triple_B': '$\\mathit{AEBP1}^{+}\\,\\mathit{COL5A2}^{+}\\,\\mathit{CA12}^{+}$'}
REPORT = []
SUFFIXES = {1: 'umap_GSE222182_primary', 2: 'umap_GSE222182_culture', 3: 'umap_GSE244889', 4: 'umap_GSE160756_tissue_background', 5: 'GSE160756_state_fraction_by_tissue', 6: 'GSE160756_frozen_signature_pseudobulk_heatmap', 7: 'umap_GSE233666_state_existence', 8: 'GSE186542_frozen_bulk_score', 9: 'state_alignment_heatmap', 10: 'independent_trajectory_scatter', 11: 'regulator_candidate_heatmap', 12: 'stratified_frozen_signature_meta_forest', 13: 'same_barcode_coexpression_effects', 14: 'ambient_proxy_and_UMI_threshold', 15: 'DIPPER_regional_protein_distribution', 16: 'leave_two_markers_out_audit', 17: 'detection_depth_standardization', 18: 'absolute_same_barcode_prevalence'}

def read(p):
    return pd.read_csv(p, sep='\t')

def cohort(s):
    return str(s).replace('GSE160756_', 'GSE160756 ').replace('_', ' ')

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

def identity_qc():
    cir = read(QC / 'leave_two_out_cohort_summary.tsv')
    sup = read(QC / 'leave_two_out_state_support.tsv')
    cir = cir[cir.endpoint.eq('lto_difference')].set_index('analysis_cohort').reindex(ORDER)
    sup = sup.set_index('analysis_cohort').reindex(ORDER)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 6.1), gridspec_kw={'width_ratios': [1.15, 1]})
    fig.subplots_adjust(left=0.18, right=0.97, bottom=0.17, top=0.8, wspace=0.34)
    ax = axes[0]
    ax.axvline(0, color='gray', ls='--', lw=0.8)
    ax.errorbar(cir.mean_difference, np.arange(8), xerr=np.vstack([cir.mean_difference - cir.ci_low, cir.ci_high - cir.mean_difference]), fmt='o', ms=4, capsize=2, color='#7B3294')
    ax.set_yticks(range(8), labels=[cohort(k) for k in ORDER])
    ax.invert_yaxis()
    ax.set_xlabel('Mean leave-two-out score difference\n(remodelling minus comparison)')
    ax.set_title('A. Score excludes ' + gene_text('AEBP1 and COL5A2'), fontsize=10.5)
    heat = sup[['other_remodeling_genes_detected_5pct', 'independent_np_genes_detected_5pct', 'score_margin_over_fixed_threshold']].astype(float)
    values = heat.copy()
    values.iloc[:, 0] /= 13
    values.iloc[:, 1] /= 8
    values.iloc[:, 2] = 0.5 + heat.iloc[:, 2] / (2 * max(heat.iloc[:, 2].abs().max(), 1e-08))
    ann = np.column_stack([heat.iloc[:, 0].map(lambda x: f'{x:.0f}/13'), heat.iloc[:, 1].map(lambda x: f'{x:.0f}/8'), heat.iloc[:, 2].map(lambda x: f'{x:.2f}')])
    sns.heatmap(values, annot=ann, fmt='', yticklabels=False, xticklabels=['Other\nremodelling\n(n/13)', 'Independent\nNP\n(n/8)', 'Score\nmargin'], cmap='YlGnBu', vmin=0, vmax=1, ax=axes[1], linewidths=0.5, cbar_kws={'label': 'Scaled support', 'shrink': 0.8})
    axes[1].set(xlabel='', ylabel='', title='B. Independent support for the frozen state')
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0, fontsize=8)
    fig.suptitle('Marker-pair exclusion sensitivity', fontsize=13, fontweight='bold', y=0.95)
    save(fig, 16)
    depth = read(QC / 'depth_standardized_cohort_summary.tsv')
    metrics = read(QC / 'depth_metrics_cohort_summary.tsv')
    metric = metrics.pivot(index='analysis_cohort', columns='metric', values='mean_log2_ratio').reindex(ORDER)
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 5.6))
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.28, top=0.79, wspace=0.47)
    axes[0].axhline(0, color='gray', ls='--', lw=0.8)
    for col, color, text in [('log2_median_umi_ratio', '#0072B2', 'Total UMI'), ('log2_median_gene_ratio', '#D55E00', 'Detected genes')]:
        axes[0].scatter(range(8), metric[col], color=color, label=text, s=28)
    axes[0].set_xticks(range(8), labels=[cohort(k) for k in ORDER], rotation=55, ha='right', fontsize=8)
    axes[0].set_ylabel('Mean log2 median ratio\n(remodelling/comparison)')
    axes[0].set_title('A. Detection-depth imbalance')
    axes[0].legend(frameon=False, fontsize=8)
    for ax, threshold in zip(axes[1:], [1, 2]):
        b = depth[depth.event.eq('double_positive') & depth.threshold_umi.eq(threshold)].set_index('analysis_cohort').reindex(ORDER)
        x = b.mean_common_support_unadjusted_risk_difference.to_numpy()
        y = b.mean_depth_standardized_risk_difference.to_numpy()
        ax.axline((0, 0), slope=1, color='gray', ls='--', lw=0.8)
        ax.scatter(x, y, s=30, color='#0072B2' if threshold == 1 else '#009E73')
        lim = max(np.nanmax(np.r_[x, y]) * 1.25, 0.4)
        ax.set(xlim=(-0.025, lim), ylim=(-0.025, lim))
        positions = [(0.37, 0.385), (0.18, 0.37), (0.055, 0.13), (0.14, 0.3), (0.18, 0.2), (0.35, 0.115), (0.32, 0.045), (0.39, 0.29)] if threshold == 1 else [(0.18, 0.34), (0.15, 0.4), (0.025, 0.075), (0.09, 0.24), (0.1, 0.13), (0.29, 0.035), (0.19, 0.01), (0.27, 0.23)]
        for i, (xx, yy) in enumerate(positions):
            ax.annotate(str(i + 1), (x[i], y[i]), xytext=(xx, yy), fontsize=9, ha='center', arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
        ax.set_xlabel('Common-support\nunadjusted risk difference')
        ax.set_ylabel('Depth-standardized risk difference')
        ax.set_title(('B' if threshold == 1 else 'C') + f'. Double positive; ' + '$\\geq$' + f' {threshold} UMI')
    fig.suptitle(gene_text('Detection-depth sensitivity of AEBP1-COL5A2 co-detection'), fontsize=13, fontweight='bold', y=0.96)
    fig.text(0.5, 0.025, '   '.join((f'{i + 1}: {cohort(k)}' for i, k in enumerate(ORDER[:4]))) + '\n' + '   '.join((f'{i + 5}: {cohort(k)}' for i, k in enumerate(ORDER[4:]))), ha='center', fontsize=8.5, linespacing=1.6)
    save(fig, 17)
    pre = read(QC / 'absolute_prevalence_cohort_summary.tsv')
    fig, axes = plt.subplots(2, 3, figsize=(12, 8.8), sharey=True)
    fig.subplots_adjust(left=0.18, right=0.98, top=0.84, bottom=0.14, hspace=0.48, wspace=0.22)
    for r, threshold in enumerate([1, 2]):
        for c, event in enumerate(EVENTS):
            ax = axes[r, c]
            b = pre[pre.event.eq(event) & pre.threshold_umi.eq(threshold)]
            for pop, off, color, lab in [('nonremodeling_NP_related', -0.12, '#999999', 'Other NP-related cells'), ('remodeling', 0.12, '#7B3294', 'Remodelling state')]:
                p = b[b.population.eq(pop)].set_index('analysis_cohort').reindex(ORDER)
                ax.errorbar(p.stratum_median_proportion, np.arange(8) + off, xerr=np.vstack([p.stratum_median_proportion - p.stratum_q1_proportion, p.stratum_q3_proportion - p.stratum_median_proportion]), fmt='o', ms=3.5, capsize=2, color=color, label=lab)
            ax.set_xscale('symlog', linthresh=0.001, linscale=0.7, base=10)
            ax.set_xlim(-0.0001, 1)
            ax.set_xticks([0, 0.001, 0.01, 0.1, 1], labels=['0', '0.1%', '1%', '10%', '100%'])
            ax.set_title(EVENT_LABELS[event] + '\n' + '$\\geq$' + f' {threshold} UMI per gene', fontsize=10)
            ax.set_xlabel('Donor/sample median prevalence')
            ax.set_yticks(range(8), labels=[cohort(k) for k in ORDER])
    axes[0, 0].invert_yaxis()
    handles, labs = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labs, loc='lower center', bbox_to_anchor=(0.5, 0.025), ncol=2, frameon=False)
    fig.suptitle('Absolute prevalence of same-barcode co-detection', fontsize=13, fontweight='bold', y=0.97)
    save(fig, 18)

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
    raw = WORK_DIR / 'Supplementary Figure S18 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S18.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        style()

        def selected_save(fig, n, note=''):
            if n == 18:
                if note:
                    fig.text(0.5, 0.012, note, ha='center', va='bottom', fontsize=8.5)
                rawsave(fig, raw, facecolor='white')
            plt.close(fig)
        save = selected_save
        identity_qc()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 18, output)
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
