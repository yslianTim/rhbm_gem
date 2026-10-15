"""Run the weak-halo diagnostic smoke check in an isolated temporary directory."""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    args = parser.parse_args()
    runner = Path(__file__).with_name('joint_offline_diagnostic.py')
    with tempfile.TemporaryDirectory(prefix='joint-offline-diagnostic-') as temporary:
        root = Path(temporary)
        subprocess.run([
            sys.executable, str(runner),
            '--work-dir', str(root / 'work'),
            '--executable', str(args.executable.resolve()),
            '--output', str(root / 'out' / 'weak-halo.json'),
        ], check=True)


if __name__ == '__main__':
    main()
