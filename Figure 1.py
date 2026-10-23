from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M1'
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
INPUT_FILES = {}
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
MAIN = PKG / 'figures/main'

def finish(fig, name, folder=MAIN):
    for ext in ['pdf', 'png']:
        fig.savefig(folder / f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def common_style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10.5, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

def figure1():
    fig, axes = plt.subplots(1, 3, figsize=(14, 6.7))
    fig.subplots_adjust(left=0.035, right=0.98, top=0.91, bottom=0.08, wspace=0.12)
    titles = ['Define the cellular expression phenotype', 'Separate the disease-associated dimensions', 'Evaluate tissue manifestation and evidence limits']
    blocks = [[('Discovery and recognition', 'GSE251686: discovery\nGSE229711 / GSE230808: independent recognition\nGSE233666 / GSE222182: tissue phenotype'), ('Identity and anatomical background', 'Fibrillar matrix plus residual cartilage expression\nLiterature FibroNP / Fibro-NPC comparison\nGSE160756: related NP / AF / CEP expression'), ('New independent supportive audit', 'GSE165722: independent clustering and fixed markers\nCandidate union and focused cluster 18 sensitivity\nConditional phenotype support; no grade-effect test')], [('Donor-level relative representation', 'GSE251686: 3 milder / 3 severe donors\nGSE153066: 8 reference / 8 degenerated donors\nGSE244889: 4 mild / 3 severe donors'), ('Expression within comparable state cells', 'Frozen 25-gene programme, unchanged directions\nRelative abundance and internal expression\nare evaluated as distinct quantities'), ('Direct donor interpretation', 'GSE244889: every donor, age context and LODO\nPositive LODO direction is not persistent significance\nGSE222182 symptoms: exploratory, not grade')], [('Four independent bulk comparisons', 'GSE56081 / GSE70362\nGSE244889 bulk / GSE186542\nNo uniform tissue-level degeneration signature'), ('Complementary and secondary layers', 'DIPPER: regional protein evidence, not colocalization\nPI3K: secondary state-internal association\nmiRNA / methylation: supplementary context'), ('Biological interpretation', 'Shared matrix features, variable cellular backgrounds\nContext-dependent relative enrichment\nNo established lineage or dominant causal mechanism')]]
    for k, (ax, title, cases) in enumerate(zip(axes, titles, blocks)):
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.text(0.01, 1.015, 'ABC'[k], fontsize=16, fontweight='bold', ha='left', va='bottom')
        ax.text(0.08, 1.015, title, fontsize=10.5, fontweight='bold', ha='left', va='bottom', wrap=True)
        for j, (head, body) in enumerate(cases):
            top = 0.955 - j * 0.315
            box = FancyBboxPatch((0.015, top - 0.265), 0.97, 0.265, boxstyle='round,pad=.008,rounding_size=.01', linewidth=0.8, edgecolor='#A8B2BB', facecolor=['#F0F5F8', '#F7F4EF', '#F1F6F3'][k])
            ax.add_patch(box)
            ax.text(0.045, top - 0.034, head, fontsize=10, fontweight='bold', va='top')
            ax.text(0.045, top - 0.095, body, fontsize=8.3, va='top', linespacing=1.65)
    fig.text(0.035, 0.025, 'Evidence roles are distinct: phenotype presence does not validate a universal disease effect or a causal pathway.', fontsize=9)
    finish(fig, 'Figure 1')

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Figure 1 raw.pdf'
    output = OUTPUT_DIR / 'Figure 1.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        common_style()
        read = cached_read
        figure1()
    finally:
        pd.read_csv, Figure.savefig = (old_read, old_save)
        plt.close('all')
    if not raw.is_file():
        raise RuntimeError('No plot output: ' + FIGURE_ID)
    output.write_bytes(raw.read_bytes())
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
