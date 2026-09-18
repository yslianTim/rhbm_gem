"""Apply the audit-only scope correction to completed, frozen paired searches.

This reproduction helper preserves the original endpoint diagnostics, then runs
the stronger ladder on every previously skipped precision target. It never runs
search, modifies fitted parameters, changes a tolerance or reuses another run's
numerical audit. New runs using the final runner perform this audit directly.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO/'tests/integration'))
import joint_abc_certification as cert


def digest(value):
    return hashlib.sha256(json.dumps(cert.coverage.joint.scientific(value), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    root, executable = (Path(p).resolve() for p in sys.argv[1:])
    cert.require(cert.read(root/'execution-status.json')['complete'], 'Original run must have completed before revision.')
    current = cert.fixed.experiment.provenance(executable, REPO)
    old = cert.read(root/'provenance.json')
    changed = [p for p, value in old['source_files'].items() if current['source_files'].get(p) != value]
    allowed = ['tests/integration/joint_abc_certification.py', 'tests/integration/joint_abc_certification_test.py',
               'tests/support/JointABCCertification.cpp']
    cert.require(sorted(changed) == sorted(allowed), 'Unexpected search-source change during audit revision.')
    revision = root.parent/(root.name+'-audit-revision')
    cert.require(not revision.exists(), 'Audit revision requires a fresh directory.')
    revision.mkdir(); selected = []
    for fit_path in sorted(root.glob('datasets/*/*/fits/*.json')):
        fit = cert.read(fit_path); audit_path = fit_path.parent.parent/'audits'/fit_path.name
        audit = cert.read(audit_path)
        if not audit.get('precision') or audit.get('ladders'): continue
        cert.require(fit['derivative_verified'], 'Unexpected skipped derivative ladder.')
        dataset = fit_path.parents[2]; copied = revision/'datasets'/dataset.name
        if not copied.exists():
            copied.mkdir(parents=True)
            for name in ('voxels.csv', 'contributors.csv', 'snapshot.json'):
                shutil.copyfile(dataset/name, copied/name)
            for variant in cert.VARIANTS:
                (copied/variant/'fits').mkdir(parents=True)
                (copied/variant/'audits').mkdir()
        target = revision/fit_path.relative_to(root); shutil.copyfile(fit_path, target)
        selected.append({'fit': str(fit_path.relative_to(root)), 'audit': str(audit_path.relative_to(root)),
                         'frozen_fit_scientific_sha256': digest(fit), 'before_audit_scientific_sha256': digest(audit)})
    cert.require(selected, 'No skipped precision targets to revise.')
    for name in ('inputs.json', 'input-hashes.json'):
        shutil.copyfile(root/name, revision/name)
    cert.audit(revision, executable)
    for row in selected:
        source = revision/row['audit']; audit = cert.read(source); fit = cert.read(root/row['fit'])
        cert.validate_audit(audit, fit)
        cert.require(len(audit['ladders']) == 3, 'Revised target is missing its full ladder.')
        cert.require(digest(fit) == row['frozen_fit_scientific_sha256'], 'Audit modified frozen fit.')
        target = root/row['audit']; previous = root/'pre-revision-audits'/Path(row['audit']).relative_to('datasets')
        previous.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(target, previous)
        shutil.copyfile(source, target)
        row['after_audit_scientific_sha256'] = digest(audit)
    shutil.copyfile(revision/'audit-provenance.json', root/'audit-provenance.json')
    cert.write(root/'audit-revision.json', {'reason': 'Require the full independent step ladder for precision targets even when the legacy two-step check passes.',
        'search_reexecuted': False, 'truth_used': False, 'changed_sources': changed, 'targets': selected,
        'seconds': cert.read(revision/'audit-execution.json')['seconds']})
    cert.write(root/'execution-status.json', {'complete': True, 'inputs_stable': True,
        'original_search_and_audit_sources_stable': True, 'revision_audit_sources_stable': True,
        'scope': 'Separate immutable search/original-audit and final targeted-audit source records.'})
    print(json.dumps({'revised_targets': len(selected), 'directory': str(root)}, indent=2))


if __name__ == '__main__': main()
