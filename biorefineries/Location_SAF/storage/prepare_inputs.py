"""Prepare storage inputs from existing, checksum-verified repository datasets."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify prepared files without copying')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    repository = root.parents[2]
    manifest = json.loads((root / 'INPUTS.json').read_text())
    errors = []
    for relative, entry in manifest.items():
        path = root / relative if args.check else repository / entry['source']
        if not path.is_file():
            errors.append(f'Missing: {path}')
        elif hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            errors.append(f'Checksum mismatch: {path}')
    if errors:
        raise SystemExit('Input verification failed:\n' + '\n'.join(errors) +
                         '\nUse the matching repository commit and include the source paths in INPUTS.json.')
    if not args.check:
        for relative, entry in manifest.items():
            destination = root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(repository / entry['source'], destination)
    print(f"{'Verified' if args.check else 'Prepared'} {len(manifest)} input files.")


if __name__ == '__main__':
    main()
