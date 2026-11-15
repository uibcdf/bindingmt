"""Capture a region, characterize its interface, and read the saved history.

Copy this directory anywhere and run with an existing compatible BindingMT
environment. The input is an authored four-ion control, not a biological case.
"""

import argparse
import hashlib
import json
from pathlib import Path

import pyunitwizard as puw

import bindingmt as bmt

FILES = ('site.json', 'region.json', 'contacts.json', 'burial.json', 'footprint.json')
EXAMPLE = 'bindingmt-site-interface-walkthrough-v1'


def _write_json(path, record):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write('\n')


def _load_records(output):
    return {
        'site': bmt.BindingSite.read(output / 'site.json', experimental=True),
        'region': bmt.BindingSite.region_from_dict(
            json.loads((output / 'region.json').read_text()), experimental=True
        ),
        'contacts': bmt.BindingInterface.read(
            output / 'contacts.json', experimental=True
        ),
        'burial': bmt.BindingInterface.read(output / 'burial.json', experimental=True),
        'footprint': bmt.ContactFootprint.read(
            output / 'footprint.json', experimental=True
        ),
    }


def _record_dicts(records):
    return {
        key: value if key == 'region' else value.to_dict()
        for key, value in records.items()
    }


def _references(records):
    return {
        key: {name: value[name] for name in ('record_id', 'content_sha256')}
        for key, value in _record_dicts(records).items()
    }


def run_example(output, *, structure_index=1):
    """Run the fixed example on one chosen frame and exclusively save its records."""
    import molsysmt as msm

    output = Path(output)
    output.mkdir()  # Refuse an existing directory before reading or calculating.
    data = Path(__file__).resolve().parent
    case = json.loads((data / 'case.json').read_text())
    source_path = data / 'source.h5msm'
    source_bytes = source_path.read_bytes()
    source_sha = hashlib.sha256(source_bytes).hexdigest()
    if source_sha != case['source']['sha256']:
        raise ValueError('The example input differs from its frozen case manifest.')
    system = msm.h5msm.read(str(source_path))
    source = {'owner': 'bindingmt.example', 'id': case['case'], 'revision': source_sha}
    recipe = case['recipe']

    site = bmt.BindingSite.from_selection(
        system,
        recipe['site_selection'],
        source=source,
        structure_indices=[structure_index],
        chemical_state='structure',
        name='Authored candidate region',
    )
    region = site.resolve_region(system, structure_index=structure_index)
    participants = {
        'region': region['projection']['exported_target_atom_indices'],
        'partner': msm.select(
            system,
            recipe['partner_selection'],
            structure_indices=[structure_index],
            chemical_state='structure',
        ),
    }
    contacts = bmt.calculate_interface_contacts(
        system,
        participants=participants,
        cutoff=recipe['contact_cutoff'],
        structure_index=structure_index,
        source_id=source['id'],
        max_dense_pairs=16,
    )
    burial = bmt.calculate_interface_burial(
        system,
        participants=participants,
        structure_index=structure_index,
        source_id=source['id'],
        probe_radius=recipe['probe_radius'],
        n_sphere_points=recipe['n_sphere_points'],
        hydrogen_policy='as_supplied',
        max_surface_evaluations=1440,
    )
    footprint = contacts.project_footprint(participant='region')
    records = dict(
        site=site, region=region, contacts=contacts, burial=burial, footprint=footprint
    )

    for key in ('site', 'contacts', 'burial', 'footprint'):
        records[key].save(output / f'{key}.json')
    _write_json(output / 'region.json', region)
    restored = _load_records(output)
    if _record_dicts(restored) != _record_dicts(records):
        raise ValueError('Saved records did not preserve the original results.')
    if source_path.read_bytes() != source_bytes:
        raise ValueError('The example input changed during the calculation.')
    with (output / 'source.h5msm').open('xb') as stream:
        stream.write(source_bytes)
    # A local example receipt links existing records; it is not a library format.
    _write_json(
        output / 'run.json',
        {
            'example': EXAMPLE,
            'source': {'file': 'source.h5msm', **source},
            'files': {
                name: hashlib.sha256((output / name).read_bytes()).hexdigest()
                for name in FILES
            },
            'records': _references(records),
        },
    )
    return restored


def read_results(output):
    """Read historical records only; never open a molecular source or recalculate."""
    output = Path(output)
    receipt = json.loads((output / 'run.json').read_text())
    if receipt['example'] != EXAMPLE:
        raise ValueError('Supply the output of this walkthrough.')
    for name in FILES:
        if (
            hashlib.sha256((output / name).read_bytes()).hexdigest()
            != receipt['files'][name]
        ):
            raise ValueError(f'Saved file differs from its run receipt: {name}')
    records = _load_records(output)
    if _references(records) != receipt['records']:
        raise ValueError('Saved record references differ from the run receipt.')
    return records


def show(records):
    """Display actual quantities explicitly in nm**2; full records keep precision."""
    axes = records['region']['projection']['axes']
    print(
        f'Structure {axes["source_structure_index"]}; chemical state {axes["chemical_state_index"]}'
    )
    print(
        'Candidate region:',
        records['region']['projection']['exported_target_atom_indices'],
    )
    print(
        'Partner:', records['contacts'].participants['partner']['source_atom_indices']
    )
    print('Contacts:', len(records['contacts'].contact_calculation['local_pairs']))
    print('Observed region footprint:', records['footprint'].source_atom_indices)
    for side, area in records['burial'].burial_areas['sides'].items():
        print(
            f'Buried area, {side}: {float(puw.get_value(area, to_unit="nm**2")):.6g} nm**2'
        )
    print(
        'Burial numerically usable:',
        records['burial'].summary()['native_surface']['burial_usable'],
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'output', type=Path, help='Fresh run directory, or saved directory with --read'
    )
    parser.add_argument(
        '--structure-index',
        type=int,
        choices=(0, 1, 2),
        default=1,
        help='Authored frame to calculate (default: 1)',
    )
    parser.add_argument(
        '--read',
        action='store_true',
        help='Read saved history without molecular input or new calculations',
    )
    args = parser.parse_args(argv)
    if args.read:
        records = read_results(args.output)
        print('Saved history; molecular source not accessed.')
    else:
        try:
            records = run_example(args.output, structure_index=args.structure_index)
        except FileExistsError:
            parser.exit(2, f'Choose a fresh output directory: {args.output}\n')
    show(records)


if __name__ == '__main__':
    main()
