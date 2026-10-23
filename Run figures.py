from pathlib import Path
import argparse
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
_bundled_python = ROOT / 'Environment/Python/python.exe'
PYTHON = _bundled_python if _bundled_python.is_file() else Path(sys.executable)
FIGURES = {**{f'M{i}': f'Figure {i}.py' for i in range(1, 8)}, **{f'S{i}': f'Supplementary Figure S{i}.py' for i in range(1, 25)}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('figures', nargs='*', default=['all'])
    parser.add_argument('--list', action='store_true')
    args = parser.parse_args()
    if args.list:
        for key, filename in FIGURES.items():
            print(f'{key}: {filename}')
        return
    requested = list(FIGURES) if args.figures == ['all'] else [key.upper() for key in args.figures]
    invalid = [key for key in requested if key not in FIGURES]
    if invalid:
        parser.error('Unknown figure: ' + ', '.join(invalid))
    environment = dict(os.environ)
    for name in ['PYTHONHOME', 'PYTHONPATH']:
        environment.pop(name, None)
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    for key in requested:
        print(key, flush=True)
        subprocess.run([str(PYTHON), '-I', '-B', '-X', 'utf8', str(ROOT / FIGURES[key])], cwd=ROOT, env=environment, check=True)


if __name__ == '__main__':
    main()
