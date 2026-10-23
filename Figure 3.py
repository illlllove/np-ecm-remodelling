from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M3'
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
INPUT_FILES = {'cep-np/08_robustness/cell_level_coexpression_audit/cohort_paired_enrichment.tsv': 'Inputs/3d47d8bd1b cohort paired enrichment.tsv', 'cep-np/08_robustness/cell_level_coexpression_audit/dipper_np_af_descriptive_contrast.tsv': 'Inputs/ae54236142 dipper np af descriptive contrast.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/depth_standardized_cohort_summary.tsv': 'Inputs/3bff1a9435 depth standardized cohort summary.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/absolute_prevalence_cohort_summary.tsv': 'Inputs/0b6331f41f absolute prevalence cohort summary.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/leave_two_out_cohort_summary.tsv': 'Inputs/1bbb6f02d1 leave two out cohort summary.tsv', 'cep-np/08_robustness/prepublication_coexpression_qc/leave_two_out_state_support.tsv': 'Inputs/dc40c112a5 leave two out state support.tsv'}
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
import shutil
import importlib.util
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
ROOT = (HERE / 'Source keys/cep-np/scripts/53_build_round2_manuscript_figures.py').resolve().parents[1]
V2 = ROOT / 'manuscript_package_v2'
MAIN = V2 / 'figures' / 'main'
BLUE = '#275D8C'
ORANGE = '#C85A0A'
GREEN = '#4E8A3A'
RED = '#B4443E'
PURPLE = '#7B3294'
GRAY = '#777777'
LIGHT_GRAY = '#B7B7B7'

def configure_style() -> None:
    sns.set_theme(style='whitegrid', context='paper')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.titlesize': 10.5, 'axes.labelsize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5, 'figure.dpi': 150, 'savefig.dpi': 400, 'axes.spines.top': False, 'axes.spines.right': False})

def label_panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.15, 1.08, label, transform=ax.transAxes, fontsize=13, fontweight='bold', va='top', ha='left')

def short_label(value: str) -> str:
    return value.replace('GSE160756_', '160756_')

def build_figure3() -> None:
    coexp = pd.read_csv(ROOT / '08_robustness' / 'cell_level_coexpression_audit' / 'cohort_paired_enrichment.tsv', sep='\t')
    depth = pd.read_csv(ROOT / '08_robustness' / 'prepublication_coexpression_qc' / 'depth_standardized_cohort_summary.tsv', sep='\t')
    prevalence = pd.read_csv(ROOT / '08_robustness' / 'prepublication_coexpression_qc' / 'absolute_prevalence_cohort_summary.tsv', sep='\t')
    lto = pd.read_csv(ROOT / '08_robustness' / 'prepublication_coexpression_qc' / 'leave_two_out_cohort_summary.tsv', sep='\t')
    support = pd.read_csv(ROOT / '08_robustness' / 'prepublication_coexpression_qc' / 'leave_two_out_state_support.tsv', sep='\t')
    dipper = pd.read_csv(ROOT / '08_robustness' / 'cell_level_coexpression_audit' / 'dipper_np_af_descriptive_contrast.tsv', sep='\t')
    order = ['GSE251686', 'GSE153066', 'GSE244889', 'GSE233666', 'GSE222182', 'GSE160756_NP', 'GSE160756_AF', 'GSE160756_CEP']
    labels = [short_label(x) for x in order]
    ypos = np.arange(len(order))[::-1]
    fig = plt.figure(figsize=(13.4, 8.7))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.08], width_ratios=[1.08, 1.0, 1.12], wspace=0.48, hspace=0.45)
    ax = fig.add_subplot(gs[0, 0])
    label_panel(ax, 'A')
    block = coexp[(coexp.event == 'double_positive') & (coexp.threshold_umi == 1)].set_index('analysis_cohort').reindex(order)
    mean = block.mean_risk_difference.to_numpy(float)
    low = block.bootstrap_mean_ci_low.to_numpy(float)
    high = block.bootstrap_mean_ci_high.to_numpy(float)
    ax.axvline(0, color=GRAY, lw=0.9, ls='--')
    ax.errorbar(mean, ypos, xerr=np.vstack([mean - low, high - mean]), fmt='o', color=BLUE, ecolor=BLUE, capsize=2.5, ms=5)
    ax.set_yticks(ypos, labels=labels)
    ax.set_xlabel('Mean risk difference\nremodelling minus comparison')
    ax.set_title('Raw-count AEBP1-COL5A2 co-detection')
    ax.grid(axis='y', visible=False)
    ax = fig.add_subplot(gs[0, 1])
    label_panel(ax, 'B')
    block = depth[(depth.event == 'double_positive') & (depth.threshold_umi == 1)].set_index('analysis_cohort').reindex(order)
    x = block.mean_common_support_unadjusted_risk_difference.to_numpy(float)
    y = block.mean_depth_standardized_risk_difference.to_numpy(float)
    flags = block.material_depth_sensitivity_flag.astype(str).str.lower().eq('true').to_numpy()
    maxv = max(np.nanmax(np.r_[x, y]) * 1.2, 0.42)
    ax.plot([0, maxv], [0, maxv], ls='--', lw=0.9, color=GRAY)
    label_positions = {'GSE153066': (0.225, 0.378, 'left'), 'GSE233666': (0.175, 0.323, 'left'), 'GSE160756_CEP': (0.365, 0.307, 'left'), 'GSE251686': (0.34, 0.267, 'left'), 'GSE222182': (0.209, 0.229, 'left'), 'GSE244889': (0.079, 0.15, 'left'), 'GSE160756_NP': (0.354, 0.116, 'left'), 'GSE160756_AF': (0.329, 0.049, 'left')}
    for idx, cohort in enumerate(order):
        edge = RED if flags[idx] else BLUE
        face = 'white' if flags[idx] else BLUE
        ax.scatter(x[idx], y[idx], s=43, facecolor=face, edgecolor=edge, linewidth=1.4, zorder=3)
        tx, ty, align = label_positions[cohort]
        ax.annotate(short_label(cohort), xy=(x[idx], y[idx]), xytext=(tx, ty), textcoords='data', ha=align, va='center', fontsize=6.4, arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.55, shrinkA=2, shrinkB=3))
    ax.set_xlim(-0.01, maxv)
    ax.set_ylim(-0.01, maxv)
    ax.set_xlabel('Common-support unadjusted difference')
    ax.set_ylabel('Depth-standardized difference')
    ax.set_title('Detection-depth sensitivity')
    ax = fig.add_subplot(gs[0, 2])
    label_panel(ax, 'C')
    block = prevalence[prevalence.population.eq('remodeling') & prevalence.threshold_umi.eq(1) & prevalence.event.isin(['triple_A', 'triple_B'])].copy()
    offsets = {'triple_A': 0.12, 'triple_B': -0.12}
    colors = {'triple_A': ORANGE, 'triple_B': PURPLE}
    names = {'triple_A': 'KRT19 triple', 'triple_B': 'CA12 triple'}
    for event in ['triple_A', 'triple_B']:
        p = block[block.event.eq(event)].set_index('analysis_cohort').reindex(order)
        med = p.stratum_median_proportion.to_numpy(float)
        lo = p.stratum_q1_proportion.to_numpy(float)
        hi = p.stratum_q3_proportion.to_numpy(float)
        yy = ypos + offsets[event]
        ax.errorbar(med, yy, xerr=np.vstack([med - lo, hi - med]), fmt='o', color=colors[event], ecolor=colors[event], capsize=2, ms=4, label=names[event])
    ax.set_xscale('symlog', linthresh=0.001, linscale=0.7, base=10)
    ax.set_xlim(-0.0001, 0.35)
    ax.set_xticks([0, 0.001, 0.01, 0.1], labels=['0', '0.1%', '1%', '10%'])
    ax.set_yticks(ypos, labels=labels)
    ax.set_xlabel('Donor/sample median prevalence')
    ax.set_title('Triple-positive cells are rare or variable')
    ax.legend(frameon=False, loc='lower right')
    ax.grid(axis='y', visible=False)
    ax = fig.add_subplot(gs[1, 0:2])
    label_panel(ax, 'D')
    block = lto[lto.endpoint.eq('lto_difference')].set_index('analysis_cohort').reindex(order)
    sup = support.set_index('analysis_cohort').reindex(order)
    mean = block.mean_difference.to_numpy(float)
    low = block.ci_low.to_numpy(float)
    high = block.ci_high.to_numpy(float)
    passed = sup.leave_two_out_support.astype(str).str.lower().eq('true').to_numpy()
    ax.axvline(0, color=GRAY, lw=0.9, ls='--')
    for idx in range(len(order)):
        face = GREEN if passed[idx] else 'white'
        edge = GREEN if passed[idx] else RED
        ax.errorbar(mean[idx], ypos[idx], xerr=[[mean[idx] - low[idx]], [high[idx] - mean[idx]]], fmt='o', mfc=face, mec=edge, mew=1.4, ecolor=edge, capsize=2.5, ms=6)
    ax.set_yticks(ypos, labels=labels)
    ax.set_xlabel('Donor/sample mean leave-two-out score difference\nremodelling minus comparison')
    ax.set_title('AEBP1 and COL5A2 are not the sole molecular support')
    ax.grid(axis='y', visible=False)
    ax.text(0.99, 0.02, 'Open symbol: complete state-support rule failed', transform=ax.transAxes, ha='right', va='bottom', fontsize=7, color=RED)
    ax = fig.add_subplot(gs[1, 2])
    label_panel(ax, 'E')
    gene_order = ['AEBP1', 'COL5A2', 'KRT19']
    ages = ['young', 'old']
    rows = []
    for gene in gene_order:
        for age in ages:
            r = dipper[(dipper.gene == gene) & (dipper.age_group == age)].iloc[0]
            rows.append((gene, age, float(r.mean_difference_NP_side_minus_AF)))
    xbase = np.arange(len(gene_order))
    width = 0.34
    for ai, age in enumerate(ages):
        vals = [next((v for g, a, v in rows if g == gene and a == age)) for gene in gene_order]
        bars = ax.bar(xbase + (ai - 0.5) * width, vals, width=width, color=BLUE if age == 'young' else LIGHT_GRAY, edgecolor='white', label=age.capitalize())
        bars[-1].set_edgecolor(RED)
        bars[-1].set_linewidth(1.3)
        bars[-1].set_hatch('//')
    ax.axhline(0, color=GRAY, lw=0.9)
    ax.set_xticks(xbase, labels=gene_order)
    ax.set_ylabel('Mean NP-side minus AF\nlog2 abundance')
    ax.set_title('DIPPER regional protein context')
    ax.legend(frameon=False, loc='upper left')
    ax.text(0.5, -0.18, 'KRT19 source protein group carries a contaminant flag', transform=ax.transAxes, ha='center', va='top', fontsize=7.2, fontweight='bold', color=RED, clip_on=False)
    fig.suptitle('Same-barcode and regional evidence support a remodelling phenotype with explicit identity limits', fontsize=13.5, fontweight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = MAIN / 'Figure3_same_barcode_identity_and_quality_audit'
    fig.savefig(out.with_suffix('.png'), bbox_inches='tight', facecolor='white')
    fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    plt.close(fig)

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
    raw = WORK_DIR / 'Figure 3 raw.pdf'
    output = OUTPUT_DIR / 'Figure 3.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        configure_style()
        build_figure3()
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
