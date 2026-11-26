"""Run a version-specific SOLMAE script with its own module search path."""
import argparse
import runpy
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('version', choices=('v4', 'v5'))
parser.add_argument('script', nargs='?')
args, extra = parser.parse_known_args()
folder = ROOT / ('solmae-fxp-integer' if args.version == 'v4' else 'solmae-full-integer')
script = args.script or ('test_fxp_solmae512.py' if args.version == 'v4' else 'test_fxp_strict_solmae512.py')
path = (folder / script).resolve()
if not path.is_file() or path.parent != folder:
    parser.error('script must exist directly in the selected version folder')
sys.path[:0] = [str(folder), str(ROOT / 'solmae-ref')]
sys.argv = [str(path), *extra]
runpy.run_path(str(path), run_name='__main__')
