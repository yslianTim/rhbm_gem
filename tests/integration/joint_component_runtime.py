#!/usr/bin/env python3
"""Fresh Map/Model runtime parity, independent geometry/replay, and repeat gates."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import os
import numpy as np
import joint_abc_components as components
import joint_abc_component_records as records
from joint_abc_profile import scientific
from observation_matching import require, read


def geometry(root):
    dataset = read(root/'dataset.json')
    dims, spacing, origin = dataset['grid_size'], dataset['spacing'], dataset['origin']
    voxels = np.atleast_2d(np.loadtxt(root/'voxels.csv', delimiter=',', skiprows=1))
    by_index = {int(row[1]): k for k, row in enumerate(voxels)}
    expected = []
    for a, atom in enumerate(dataset['atoms']):
        p = atom['position']
        lo = [max(0, min(dims[k], math.floor((p[k]-2.5-origin[k])/spacing[k]))) for k in range(3)]
        hi = [max(-1, min(dims[k]-1, math.floor((p[k]+2.5-origin[k])/spacing[k]))) for k in range(3)]
        for z in range(lo[2], hi[2]+1):
            for y in range(lo[1], hi[1]+1):
                for x in range(lo[0], hi[0]+1):
                    position = [math.fma(float(i), h, o) for i, h, o in zip((x, y, z), spacing, origin)]
                    dx, dy, dz = [v-c for v, c in zip(position, p)]
                    square = math.fma(dz, dz, math.fma(dy, dy, dx*dx))
                    if square <= 6.25:
                        index = x+dims[0]*(y+dims[1]*z)
                        require(index in by_index, 'Runtime ROI lost a contributing voxel.')
                        row = by_index[index]
                        require(voxels[row, 2:5].tolist() == position, 'Runtime coordinate arithmetic differs.')
                        expected.append((row, a, square))
    data = records.load(root)
    require(sorted(expected) == list(map(tuple, data['table'])), 'Runtime support differs from independent geometry.')
    require(set(row for row, _, _ in expected) == set(range(len(voxels))), 'Runtime union contains an extraneous row.')
    records.verify_census(data, read(root/'census.json'))
    return data


def summarize(root):
    completion = read(root/'completion.json')
    roots = [root/name for name in completion['datasets']] if 'datasets' in completion else [root]
    pairs, local_records = [], []
    for dataset in roots:
        data = geometry(dataset); census = read(dataset/'census.json')
        done = read(dataset/'completion.json')
        for case in done['cases']:
            target = dataset/'cases'/case
            y = data['y64' if case.endswith('double') else 'y32']
            fits = [read(target/(name+'-fit.json')) for name in ('monolithic', 'assembled')]
            certificates = []
            for name in ('monolithic', 'assembled'):
                audit_dir = target/name
                fit, audit = read(audit_dir/'fit.json'), read(audit_dir/'audit.json')
                components.validate_endpoint_scope(fit, read(audit_dir/'context.json'))
                cert = components.certify(data, y, fit, audit, 'assembled-global' if name == 'assembled' else 'monolithic-global')
                components.write(audit_dir/'certificate.json', cert); certificates.append(cert)
            pair = components.endpoint_parity(data, y, *fits, certificates, done['requires_regular_parity'])
            require(pair['passed'], f'Fresh endpoint parity failed: {dataset.name}/{case}')
            runtime = read(target/'runtime.json')
            require(runtime['regular_certificate'] == 'not-run', 'Runtime claimed an offline certificate.')
            require(runtime['prediction_reassembly_difference'] == 0, 'Runtime coefficients changed during reassembly.')
            pairs.append(dict(dataset=dataset.name, case=case, **pair))
            for k, component in enumerate(census['components']):
                local = components.subset(data, component); local_y = local['y64' if case.endswith('double') else 'y32']
                directory = target/'local-components'/str(k)
                fit, audit, scope = (read(directory/name) for name in ('fit.json', 'audit.json', 'scope.json'))
                if scope['weak_direction_available']: components.validate_endpoint_scope(fit, read(directory/'context.json'))
                cert = components.certify(local, local_y, fit, audit, 'component-local')
                cert['checks'].update(local_weak_direction=scope['weak_direction_available'],
                    actual_state_profile=fit.get('assembled_profile_agrees', False), usable_state=fit['usable_state'])
                cert['failures'] = [k for k, v in cert['checks'].items() if v is False]
                cert['unavailable'] = [k for k, v in cert['checks'].items() if v is None]
                cert['regular_qualified'] = all(v is True for v in cert['checks'].values())
                components.write(directory/'certificate.json', cert)
                local_records.append(dict(dataset=dataset.name, case=case, component=component['id'], regular=cert['regular_qualified']))
    result = dict(passed=True, paired_branches=len(pairs), regular_pairs=sum(p['both_regular'] for p in pairs),
                  local_certificates=len(local_records), local_regular=sum(p['regular'] for p in local_records),
                  geometry_verified=True, historical_files_required=False)
    components.write(root/'summary.json', result)
    components.csv_write(root/'endpoint-parity.csv', pairs)
    components.csv_write(root/'local-certificates.csv', local_records)
    return result


def run(args):
    output, executable = args.output.resolve(), args.executable.resolve()
    require(not output.exists(), 'Use a fresh runtime output directory.')
    before = components.provenance(executable)
    command = [str(executable), 'joint-component-physical', str(output)] if args.command == 'physical' else [
        str(executable), 'joint-component-runtime', str(args.model.resolve()), str(args.map.resolve()), str(output)]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_suffix('.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, env=dict(os.environ, **components.NUMERIC_ENV))
    require(before == components.provenance(executable), 'Runtime source or executable changed during execution.')
    components.write(output/'provenance.json', before)
    print(json.dumps(summarize(output)))


def compare(left, right):
    excluded = {'provenance.json', 'comparison.json'}
    a = {p.relative_to(left) for p in left.rglob('*.json') if p.name not in excluded}
    b = {p.relative_to(right) for p in right.rglob('*.json') if p.name not in excluded}
    differences = [str(p) for p in a ^ b]
    differences += [str(p) for p in sorted(a & b) if scientific(read(left/p)) != scientific(read(right/p))]
    require(not differences, 'Runtime scientific repeat differs: '+str(differences))
    return dict(passed=True, scientific_records=len(a), differences=[])


def frozen_parity(reference, run, reference_audits):
    differences, counts = [], {}
    datasets = set(components.certification.coverage.DATASETS) | set(components.COMPOSITES)
    for scope, old, new in (("search", reference/"datasets", run/"datasets"),
                            ("audit", reference_audits, run/"audits")):
        paths = [p for p in old.rglob("*.json") if p.name != "completion.json" and
                 p.relative_to(old).parts[0] in datasets]
        require(bool(paths), "Missing frozen "+scope+" records.")
        counts[scope] = len(paths)
        for path in paths:
            relative = path.relative_to(old); target = new/relative
            if not target.exists() or scientific(read(path)) != scientific(read(target)):
                differences.append(scope+"/"+str(relative))
    return dict(passed=not differences, compared_records=counts, differences=differences,
                reference="pre-extraction exact-component records",
                excluded=["elapsed seconds", "process peak RSS", "process completion manifests"],
                additional_component_local_records="separate; never substituted for historical scopes")


def main():
    parser = argparse.ArgumentParser(description=__doc__); commands = parser.add_subparsers(dest='command', required=True)
    for name in ('run', 'physical'):
        sub = commands.add_parser(name)
        sub.add_argument('--executable', type=Path, required=True); sub.add_argument('--output', type=Path, required=True)
        if name == 'run':
            sub.add_argument('--model', type=Path, required=True); sub.add_argument('--map', type=Path, required=True)
    sub = commands.add_parser('summarize'); sub.add_argument('--run', type=Path, required=True)
    sub = commands.add_parser('compare')
    for name in ('left', 'right', 'output'): sub.add_argument('--'+name, type=Path, required=True)
    sub = commands.add_parser("frozen-parity")
    for name in ("reference", "run", "reference-audits", "output"): sub.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command in ('run', 'physical'): run(args)
    elif args.command == 'summarize': print(json.dumps(summarize(args.run)))
    elif args.command == "frozen-parity":
        result = frozen_parity(args.reference, args.run, args.reference_audits)
        components.write(args.output, result); print(json.dumps(result))
        require(result["passed"], "Frozen scientific records differ.")
    else:
        result = compare(args.left, args.right); components.write(args.output, result); print(json.dumps(result))


if __name__ == '__main__': main()
