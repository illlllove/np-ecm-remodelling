from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'S20'
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
INPUT_FILES = {'cep-np/08_robustness/literature_phenotype_context_audit/clinical_context_descriptive_profiles.tsv': 'Inputs/507e457307 clinical context descriptive profiles.tsv', 'cep-np/08_robustness/literature_phenotype_context_audit/programme_paired_descriptive_contrasts.tsv': 'Inputs/f08f40adc4 programme paired descriptive contrasts.tsv', 'cep-np/08_robustness/literature_phenotype_context_audit/cohort_codetection_median_iqr.tsv': 'Inputs/d03832b02b cohort codetection median iqr.tsv', 'cep-np/08_robustness/literature_phenotype_context_audit/cohort_gene_expression_medians.tsv': 'Inputs/6948a5e389 cohort gene expression medians.tsv'}
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
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import seaborn as sns
ROOT = (HERE / 'Source keys/cep-np/scripts/62_literature_phenotype_context_audit.py').resolve().parents[1]
OUT = ROOT / '08_robustness/literature_phenotype_context_audit'
PROGRAMS = {'Fibrillar matrix': ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2'], 'Matrix processing': ['POSTN', 'PCOLCE', 'MMP2', 'MMP14'], 'Cartilage identity': ['ACAN', 'COL2A1', 'SOX9', 'COMP', 'CILP'], 'NP-associated identity': ['KRT8', 'KRT18', 'KRT19', 'CA12'], 'Published FibroNP context': ['COL1A1', 'COL3A1', 'POSTN', 'ASPN'], 'Published Fibro-NPC matrix context': ['COL1A1', 'MMP2', 'MMP9', 'MMP13', 'ADAMTS2', 'ADAMTS4'], 'Published FAP-related context': ['FAP', 'S100A4', 'COL1A1'], 'Immune context': ['PTPRC', 'ITGAM', 'LST1', 'TYROBP', 'FCER1G', 'CD68', 'LYZ'], 'Contractile context': ['ACTA2', 'TAGLN', 'MYL9']}
POPS = ['Frozen remodelling', 'Other NP/mesenchymal', 'Immune-positive control', 'Other excluded']
CONTEXTS = ['GSE244889', 'GSE233666', 'GSE222182', 'GSE160756 NP', 'GSE160756 AF', 'GSE160756 CEP']
MIN_CELLS = 20

def savefig(fig, name):
    for suffix in ['.png', '.pdf']:
        fig.savefig(OUT / (name + suffix), dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)

def plots(cl, pair, es, gs):
    sns.set_theme(style='ticks', context='paper')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'pdf.fonttype': 42, 'mathtext.fontset': 'dejavusans'})
    fig = plt.figure(figsize=(14.5, 12))
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 1.8], left=0.17, right=0.94, top=0.91, bottom=0.14, hspace=0.55)
    ax = fig.add_subplot(grid[0])
    programs = list(PROGRAMS)
    mat = pair.pivot(index='context', columns='program', values='mean_state_minus_other_NP').reindex(index=CONTEXTS, columns=programs)
    sns.heatmap(mat, annot=True, fmt='.2f', center=0, cmap='vlag', ax=ax, cbar_kws={'label': 'Mean paired difference\nwithin-study descriptive z coordinate'})
    ax.set_xticklabels([s.replace('Published ', '').replace(' context', '').replace(' identity', '\nidentity').replace('Fibro-NPC matrix', 'Fibro-NPC\nmatrix') for s in programs], rotation=25, ha='right')
    ax.set(xlabel='', ylabel='', title='A  Frozen remodelling state versus other NP/mesenchymal cells')
    ax = fig.add_subplot(grid[1])
    plotted = ['COL1A1', 'COL3A1', 'COL5A2', 'POSTN', 'PCOLCE', 'MMP2', 'ACAN', 'COL2A1', 'SOX9', 'COMP', 'KRT19', 'CA12', 'SFRP4', 'FAP', 'S100A4', 'ASPN', 'FBLN1', 'PTPRC', 'ITGAM', 'LST1', 'TYROBP', 'ACTA2', 'TAGLN', 'AEBP1']
    rowkeys = [(c, p) for c in CONTEXTS for p in POPS[:3]]
    mx = gs.set_index(['context', 'population', 'gene'])
    norm = Normalize(0, 4.0)
    cmap = plt.get_cmap('viridis')
    for y, (context, pop) in enumerate(rowkeys):
        for x, g in enumerate(plotted):
            if (context, pop, g) not in mx.index:
                continue
            r = mx.loc[context, pop, g]
            ax.scatter(x, y, s=95 * r.median_detected_fraction, color=cmap(norm(r.median_mean_log1p_cp10k)), edgecolor='none')
    ax.set_xlim(-0.6, len(plotted) - 0.4)
    ax.set_ylim(len(rowkeys) - 0.4, -0.6)
    ax.set_xticks(range(len(plotted)), labels=plotted, rotation=65, ha='right', fontstyle='italic')
    ax.set_yticks(range(len(rowkeys)), labels=[f'{c} | ' + {'Frozen remodelling': 'State', 'Other NP/mesenchymal': 'Other NP', 'Immune-positive control': 'Immune control'}[p] for c, p in rowkeys])
    for y in np.arange(2.5, len(rowkeys), 3):
        ax.axhline(y, color='#DDDDDD', lw=0.7)
    ax.set_title('B  Donor/library-balanced gene expression and detection', loc='left', pad=13)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.018, label='Median mean log1p CP10k')
    handles = [ax.scatter([], [], s=95 * v, color='#555555', label=f'{100 * v:.0f}%') for v in [0.1, 0.5, 1]]
    fig.legend(handles=handles, title='Median detection', loc='lower center', bbox_to_anchor=(0.55, 0.045), ncol=3, frameon=False)
    sns.despine(ax=ax)
    fig.suptitle('Shared matrix programmes and distinct expression backgrounds', fontsize=15, fontweight='bold', y=0.98)
    fig.text(0.5, 0.006, 'Marker context is descriptive, not subtype transfer or independent cell-lineage validation. Blank rows: no evaluable positive-control stratum.', ha='center', fontsize=8.5)
    savefig(fig, 'Figure_L1_literature_marker_context')
    fig = plt.figure(figsize=(13.5, 10.7))
    grid = fig.add_gridspec(2, 6, left=0.14, right=0.98, top=0.9, bottom=0.17, hspace=0.75, wspace=0.7)
    events = ['PTPRC_COL1A1', 'FAP_COL1A1', 'SFRP4_COL1A1', 'COL1A1_cartilage_ge1', 'COL1A1_cartilage_ge2']
    labs = ['$\\mathit{PTPRC}$ + $\\mathit{COL1A1}$', '$\\mathit{FAP}$ + $\\mathit{COL1A1}$', '$\\mathit{SFRP4}$ + $\\mathit{COL1A1}$', '$\\mathit{COL1A1}$ + cartilage ($\\geq$1/3)', '$\\mathit{COL1A1}$ + cartilage ($\\geq$2/3)']
    for k, threshold in enumerate([1, 2]):
        ax = fig.add_subplot(grid[0, k * 3:k * 3 + 3])
        b = es[es.population.eq(POPS[0]) & es.threshold_umi.eq(threshold)]
        m = b.pivot(index='context', columns='event', values='median_fraction').reindex(index=CONTEXTS, columns=events) * 100
        sns.heatmap(m, annot=True, fmt='.2f', cmap='YlGnBu', vmin=0, vmax=100, ax=ax, cbar_kws={'label': 'Median positive cells (%)', 'shrink': 0.8})
        ax.set_xticklabels(labs, rotation=35, ha='right', fontsize=8)
        ax.set(xlabel='', ylabel='', title=f'{chr(65 + k)}  Raw co-detection in the frozen state; ' + '$\\geq$' + f' {threshold} UMI/gene')
    palette = {'MDD': '#4C78A8', 'SDD': '#E1812C', 'asymptomatic': '#4C78A8', 'back_pain': '#E1812C', '': '#777777'}
    for k, study in enumerate(['GSE244889', 'GSE233666', 'GSE222182']):
        ax = fig.add_subplot(grid[1, k * 2:k * 2 + 2])
        b = cl[cl.study.eq(study) & (cl.state_cells >= MIN_CELLS)]
        ycol = 'Fibrillar matrix_mean_log1p_cp10k'
        for group, d in b.groupby('condition', dropna=False):
            group = str(group)
            name = {'MDD': 'Mild', 'SDD': 'Severe', 'asymptomatic': 'Asymptomatic', 'back_pain': 'Back pain', '': 'Herniated NP'}.get(group, group)
            ax.scatter(d.state_fraction, d[ycol], label=name, color=palette.get(group, '#777777'), s=35)
        ax.set_title(f'{chr(67 + k)}  {study}')
        ax.set_xlabel('Frozen state fraction\nwithin broad NP/mesenchymal cells')
        ax.set_ylabel('Within-state fibrillar-matrix\nmean log1p CP10k')
        ax.legend(frameon=False, fontsize=7.5, loc='best')
        sns.despine(ax=ax)
    fig.suptitle('Mixed phenotype prevalence and clinical-context description', fontsize=15, fontweight='bold', y=0.98)
    fig.text(0.5, 0.035, 'Cartilage markers: $\\mathit{ACAN}$, $\\mathit{COL2A1}$, $\\mathit{SOX9}$. Every required gene must meet the stated UMI threshold.' + '\nBottom panels are descriptive, unadjusted donor plots. Their expression coordinate is not the frozen 25-gene score.', ha='center', fontsize=9, linespacing=1.5)
    savefig(fig, 'Figure_L2_codetection_and_clinical_context')

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
    raw = WORK_DIR / 'Supplementary Figure S20 raw.pdf'
    output = OUTPUT_DIR / 'Supplementary Figure S20.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:

        def selected_save(fig, name):
            if name == 'Figure_L2_codetection_and_clinical_context':
                rawsave(fig, raw, dpi=300, facecolor='white', bbox_inches='tight')
            plt.close(fig)
        savefig = selected_save
        root = Path(MANIFEST['canonicalRoot']).parent / '08_robustness/literature_phenotype_context_audit'
        tables = [cached_read(root / name, sep='\t') for name in ['clinical_context_descriptive_profiles.tsv', 'programme_paired_descriptive_contrasts.tsv', 'cohort_codetection_median_iqr.tsv', 'cohort_gene_expression_medians.tsv']]
        tables[0]['condition'] = tables[0]['condition'].fillna('')
        plots(*tables)
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    wrap_supplement(raw, 20, output)
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
