from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M7'
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
INPUT_FILES = {'cep-np/09_directiona_composition_bulk/07_regulatory/state_internal_conditional_associations.tsv': 'Inputs/0dbd9c467e state internal conditional associations.tsv'}
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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
ROOT = (HERE / 'Source keys/cep-np/scripts/38_directionA_g6_manuscript_integration.py').resolve().parents[1]
BASE = ROOT / '09_directionA_composition_bulk'
FIG = BASE / '08_figures/final_main'

def figure6():
    a = pd.read_csv(BASE / '07_regulatory/state_internal_conditional_associations.tsv', sep='\t')
    a = a.loc[~a.predictor_type.eq('AEBP1_expression_self_included_sensitivity')].copy()
    a['label'] = a.predictor.str.replace('_activity', '', regex=False).str.replace('_log1p_CPM', '', regex=False)
    order = ['PI3K', 'TGFb', 'SMAD3', 'AEBP1']
    models = ['unadjusted_spearman', 'partial_spearman_control_condition', 'partial_spearman_control_remodeling_fraction']
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for ax, study in zip(axes, ['GSE251686', 'GSE153066']):
        q = a.loc[a.study_id.eq(study)].pivot(index='label', columns='model', values='rho').reindex(index=order, columns=models)
        sns.heatmap(q, vmin=-1, vmax=1, cmap='vlag', center=0, annot=True, fmt='.2f', cbar=ax is axes[-1], ax=ax, xticklabels=['Unadjusted', 'Control condition', 'Control composition'])
        ax.set(title=study, xlabel='', ylabel='')
        ax.tick_params(axis='x', rotation=30)
        ax.tick_params(axis='y', rotation=0)
    fig.tight_layout()
    fig.savefig(FIG / 'Figure6_state_internal_regulatory_associations.png', dpi=300)
    fig.savefig(FIG / 'Figure6_state_internal_regulatory_associations.pdf')

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Figure 7 raw.pdf'
    output = OUTPUT_DIR / 'Figure 7.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        sns.set_theme(style='whitegrid', context='notebook')
        figure6()
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
