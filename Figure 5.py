from __future__ import annotations
from pathlib import Path
import io, json, math, os, re, shutil, subprocess, sys
HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUTPUT_DIR = HERE / 'Outputs'
FIGURE_ID = 'M5'
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
INPUT_FILES = {'cep-np/09_directiona_composition_bulk/05_counterfactual_bulk/counterfactual_grid.tsv': 'Inputs/a290f5e897 counterfactual grid.tsv', 'cep-np/09_directiona_composition_bulk/06_deconvolution/deconvolution_benchmark.tsv': 'Inputs/2e113c9ea4 deconvolution benchmark.tsv', 'cep-np/09_directiona_composition_bulk/06_deconvolution/deconvolution_leave_one_donor_metrics.tsv': 'Inputs/bd9c6741c0 deconvolution leave one donor metrics.tsv'}
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
COL = {'GSE251686': '#1f4e79', 'GSE153066': '#c55a11'}

def figure5():
    grid = pd.read_csv(BASE / '05_counterfactual_bulk/counterfactual_grid.tsv', sep='\t')
    bench = pd.read_csv(BASE / '06_deconvolution/deconvolution_benchmark.tsv', sep='\t').query("scope=='observed_range'")
    donor = pd.read_csv(BASE / '06_deconvolution/deconvolution_leave_one_donor_metrics.tsv', sep='\t')
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, (study, q) in zip(axes[0], grid.groupby('study_id', sort=False)):
        ax.plot(q.disease_counterfactual_fraction, q.disease_minus_reference_score, color=COL[study], lw=2)
        ax.fill_between(q.disease_counterfactual_fraction, q.bootstrap_ci_low, q.bootstrap_ci_high, color=COL[study], alpha=0.18)
        ax.axhline(0, color='black', lw=1)
        ax.axvline(q.disease_observed_mean_fraction.iloc[0], color='#c55a11', ls='--')
        ax.axvline(q.reference_observed_mean_fraction.iloc[0], color='#548235', ls=':')
        ax.set(title=f'{study} observed-range counterfactual', xlabel='Remodeling-state fraction', ylabel='Disease minus reference score')
    x = np.arange(len(bench))
    axes[1, 0].bar(x - 0.18, bench.mae, 0.36, color='#5b9bd5', label='MAE')
    axes[1, 0].bar(x + 0.18, bench.spearman_rho, 0.36, color='#ed7d31', label='Spearman rho')
    axes[1, 0].axhline(0.1, color='#5b9bd5', ls='--', lw=1)
    axes[1, 0].axhline(0.6, color='#ed7d31', ls='--', lw=1)
    axes[1, 0].set_xticks(x, bench.study_id)
    axes[1, 0].set(title='Prespecified deconvolution gate', ylabel='Metric')
    axes[1, 0].legend(frameon=False)
    sns.boxplot(data=donor, x='study_id', y='mae', hue='study_id', palette=COL, legend=False, ax=axes[1, 1], showfliers=False)
    sns.stripplot(data=donor, x='study_id', y='mae', color='black', size=5, jitter=0.12, ax=axes[1, 1])
    axes[1, 1].axhline(0.1, color='black', ls='--', lw=1)
    axes[1, 1].set(title='Held-out donor calibration error', ylabel='Donor-specific MAE', xlabel='')
    for i, ax in enumerate(axes.ravel()):
        ax.text(-0.12, 1.04, chr(65 + i), transform=ax.transAxes, fontweight='bold', fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG / 'Figure5_counterfactual_and_deconvolution_gate.png', dpi=300)
    fig.savefig(FIG / 'Figure5_counterfactual_and_deconvolution_gate.pdf')

def rawsave(fig, target, **kwargs):
    kwargs.pop('format', None)
    ORIGINALSAVE(fig, target, **kwargs)

def make_figure():
    global read, save, picture, UMAPS, savefig, hedges_g_ci
    np.random.seed(20261020)
    plt.rcdefaults()
    plt.close('all')
    raw = WORK_DIR / 'Figure 5 raw.pdf'
    output = OUTPUT_DIR / 'Figure 5.pdf'
    old_read, old_save = (pd.read_csv, Figure.savefig)

    def onlysave(fig, path, *args, **kwargs):
        if Path(path).suffix.lower() == '.pdf':
            rawsave(fig, raw, **kwargs)
    pd.read_csv = cached_read
    Figure.savefig = onlysave
    try:
        sns.set_theme(style='whitegrid', context='notebook')
        figure5()
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
