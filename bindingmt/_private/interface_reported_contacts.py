"""Owned reduction of admitted original PISA endpoints; no molecular kernels."""

import math
from collections.abc import Mapping

import numpy as np
import pyunitwizard as puw
from pyunitwizard.record import QuantityRecord

from bindingmt._private.interface_records import _require
from bindingmt._private.interface_topomt import (
    GROUPS,
    _fields,
    _sha,
    _text,
    validate_observation,
    validate_receiver,
)
from bindingmt._private.smonitor import CapabilityMismatchError

DEFINITION = 'bindingmt.pisa_reported_contact_endpoints@1'
SUPPORTED_GROUPS = ('h-bonds', 'other-bonds')
MAX_SOURCE_ATOMS = 4096
CHECKS = [
    'bounded complete original native MolSys',
    'original PISA artifact/run/context recovered twice through TopoMT',
    'complete original participant chains queried through MolSysMT',
    'whole source frame geometry unchanged before/after',
]
PRECONDITIONS = [
    'complete intended original molecular identity and reference state',
    'borrowed source identity/chemistry/state and geometry unchanged throughout the call',
    'no concurrent source or artifact mutation; serial witnesses are not atomicity',
]


def selected_groups(groups):
    _require(
        isinstance(groups, (list, tuple))
        and bool(groups)
        and all(_text(g) for g in groups)
        and len(set(groups)) == len(groups),
        'contact_groups',
        'Select a nonempty ordered list of distinct original contact groups.',
    )
    _require(
        all(g in GROUPS for g in groups),
        'contact_groups',
        'Unknown original PISA contact group.',
    )
    unsupported = [g for g in groups if g not in SUPPORTED_GROUPS]
    if unsupported:
        raise CapabilityMismatchError(
            capability=f'PISA groups {unsupported}',
            engine='bindingmt.reported_contacts',
        )
    return list(groups)


def finding(observation, identity):
    _require(
        observation['outcome'] == 'completed' and _text(identity),
        'finding_id',
        'Select an original finding of a completed analysis.',
    )
    matches = [f for f in observation['findings'] if f['source_id'] == identity]
    _require(
        len(matches) == 1,
        'finding_id',
        'The completed analysis has no such interface finding.',
    )
    return matches[0]


def aggregate_endpoints(groups, source_to_local):
    """Union already mapped observations, preserving every group/row/endpoint witness.

    This is a binding-evidence reduction, not a molecular correspondence solver.
    The caller supplies the admitted map and retains its own authority/context.
    """
    _require(
        isinstance(source_to_local, Mapping)
        and all(
            type(s) is int and s >= 0 and type(t) is int and t >= 0
            for s, t in source_to_local.items()
        )
        and len(set(source_to_local.values())) == len(source_to_local),
        'mapping',
        'Use explicit unique nonnegative source/local labels.',
    )
    witnesses = {}
    for group, rows in groups:
        for row in rows:
            _fields(
                row,
                {'source_atom_index', 'contact_index', 'endpoint_index'},
                'contributors',
            )
            source = row['source_atom_index']
            _require(
                type(source) is int
                and source in source_to_local
                and type(row['contact_index']) is int
                and row['contact_index'] >= 0
                and row['endpoint_index'] in (0, 1)
                and type(row['endpoint_index']) is int,
                'contributors',
                'Retain explicit mapped source endpoints and row positions.',
            )
            witnesses.setdefault(source, []).append(
                {
                    'group': group,
                    'contact_index': row['contact_index'],
                    'endpoint_index': row['endpoint_index'],
                }
            )
    members = sorted(witnesses)
    return {
        'source_atom_indices': members,
        'local_atom_indices': [source_to_local[s] for s in members],
        'contributors': [witnesses[s] for s in members],
        'value_state': 'evaluated_nonempty' if members else 'evaluated_empty',
    }


def distance_field(group):
    return f'bindingmt.topomt.pisa.reported_distances.{group}'


def distance_values(rows):
    values = []
    for row in rows:
        packet = row['distance']
        value = packet.get('value')
        _require(
            packet.get('unit') == 'angstroms'
            and type(value) in (int, float)
            and math.isfinite(value)
            and value >= 0,
            'distances',
            'Preserve finite original reported distances and their explicit angstrom unit.',
        )
        values.append(value)
    return np.asarray(values, dtype=np.float64)


def distance_records(finding_record, groups):
    return {
        group: QuantityRecord.from_quantity(
            puw.quantity(
                distance_values(finding_record['fields']['contacts'][group]),
                'angstrom',
                form='pint',
            ),
            field=distance_field(group),
        ).to_dict()
        for group in groups
    }


def payload(projection):
    original = finding(projection['observation'], projection['finding_id'])
    chains = [p['chain_id'] for p in original['fields']['participants']]
    side = chains.index(
        projection['participant_bindings'][projection['participant']]['owner_chain_id']
    )
    axis = set(
        projection['participant_bindings'][projection['participant']][
            'source_atom_indices'
        ]
    )
    groups = []
    coverage = []
    for group in projection['contact_groups']:
        rows = original['fields']['contacts'][group]
        endpoints = []
        for index, row in enumerate(rows):
            pair = row.get('source_atom_indices')
            _require(
                isinstance(pair, list)
                and len(pair) == 2
                and pair[side] in axis
                and row['original'].get(f'chain-{side + 1}') == chains[side],
                'mapping',
                'Original endpoint/participant correspondence is unavailable or contradictory.',
            )
            endpoints.append(
                {
                    'source_atom_index': pair[side],
                    'contact_index': index,
                    'endpoint_index': side,
                }
            )
        groups.append((group, endpoints))
        coverage.append(
            {
                'name': group,
                'row_count': len(rows),
                'value_state': 'evaluated_nonempty' if rows else 'evaluated_empty',
            }
        )
    return {
        'participant': projection['participant'],
        'source_id': projection['source_id'],
        'contact_groups': coverage,
        **aggregate_endpoints(groups, dict(projection['local_atom_map'])),
    }


def validate_projection(value):
    _fields(
        value,
        {
            'definition',
            'analysis_name',
            'analysis_ref',
            'parent_ref',
            'observation',
            'finding_id',
            'source_id',
            'participant',
            'participant_bindings',
            'contact_groups',
            'local_atom_map',
            'source_context',
            'reported_distances',
            'receiver',
        },
        'reported_projection',
    )
    _require(
        value['definition'] == DEFINITION
        and _text(value['analysis_name'])
        and _text(value['source_id']),
        'reported_projection',
        'Preserve the admitted reported-contact definition and identities.',
    )
    validate_receiver(value['receiver'], 'BindingInterface.project_reported_footprint')
    validate_observation(value['observation'])
    observation = value['observation']
    original = finding(observation, value['finding_id'])
    _require(
        value['analysis_ref']
        == {k: observation[k] for k in ('provider', 'backend', 'run_id')},
        'reported_projection',
        'Preserve the exact analysis/run reference.',
    )
    _fields(value['parent_ref'], {'record_id', 'content_sha256'}, 'parent_ref')
    _require(
        _text(value['parent_ref']['record_id'])
        and _sha(value['parent_ref']['content_sha256']),
        'parent_ref',
        'Preserve the exact original parent identity/digest.',
    )
    _require(
        selected_groups(value['contact_groups']) == value['contact_groups'],
        'contact_groups',
        'Preserve the original request order.',
    )
    context = value['source_context']
    state_context = observation['metadata'].get('chemical_state_context')
    _require(
        isinstance(state_context, Mapping),
        'source_context',
        'A historical molecular check requires known original state correspondence.',
    )
    _fields(
        context,
        {
            'scope',
            'context_policy',
            'source_atom_count',
            'source_structure_index',
            'chemical_state_index',
            'geometry_sha256_nm_float64',
            'bundle_sha256',
            'performed_checks',
            'caller_preconditions',
            'assessment_kind',
        },
        'source_context',
    )
    count = context['source_atom_count']
    _require(
        context['scope'] == 'molecular_correspondence'
        and context['context_policy'] == 'require_match'
        and type(count) is int
        and 0 < count <= MAX_SOURCE_ATOMS
        and context['assessment_kind'] == 'historical_check_receipt',
        'source_context',
        'Preserve the scoped historical source-correspondence check, not a current finding guarantee.',
    )
    _require(
        context['source_structure_index']
        == observation['metadata'].get('structure_indices')
        and type(context['source_structure_index']) is int
        and context['source_structure_index'] >= 0
        and context['chemical_state_index'] == state_context.get('index')
        and type(context['chemical_state_index']) is int
        and context['chemical_state_index'] >= 0
        and all(i < count for i in observation['metadata']['selected_atom_indices'])
        and _sha(context['geometry_sha256_nm_float64'])
        and _sha(context['bundle_sha256']),
        'source_context',
        'Preserve exact original frame/state and checked artifact/geometry identities.',
    )
    _require(
        context['performed_checks'] == CHECKS
        and context['caller_preconditions'] == PRECONDITIONS,
        'source_context',
        'Preserve performed facts separately from unchanged-input caller preconditions.',
    )
    _require(
        isinstance(value['local_atom_map'], list)
        and all(
            isinstance(pair, list)
            and len(pair) == 2
            and all(type(i) is int for i in pair)
            for pair in value['local_atom_map']
        )
        and value['local_atom_map'] == [[i, i] for i in range(count)],
        'mapping',
        'This first molecular profile uses the complete original identity axis.',
    )
    bindings = value['participant_bindings']
    _require(
        isinstance(bindings, Mapping)
        and len(bindings) == 2
        and _text(value['participant'])
        and value['participant'] in bindings,
        'participant',
        'Select one of the two admitted original participant bindings.',
    )
    labels, axes = [], []
    for name, binding in bindings.items():
        _fields(binding, {'owner_chain_id', 'source_atom_indices'}, 'participants')
        indices = binding['source_atom_indices']
        _require(
            _text(name)
            and _text(binding['owner_chain_id'])
            and isinstance(indices, list)
            and indices
            and all(type(i) is int and 0 <= i < count for i in indices)
            and indices == sorted(set(indices)),
            'participants',
            'Preserve complete original chain membership and exact source indices.',
        )
        labels.append(binding['owner_chain_id'])
        axes.append(set(indices))
    _require(
        set(labels) == {p['chain_id'] for p in original['fields']['participants']}
        and not axes[0] & axes[1]
        and original['atom_indices'] is not None,
        'mapping',
        'Require two distinct original chains and available original endpoint maps, even for selected empty groups.',
    )
    _fields(value['reported_distances'], value['contact_groups'], 'distances')
    chain_axes = {
        b['owner_chain_id']: set(b['source_atom_indices']) for b in bindings.values()
    }
    chains = [p['chain_id'] for p in original['fields']['participants']]
    for group in value['contact_groups']:
        _require(
            all(
                all(
                    row['source_atom_indices'][side] in chain_axes[chain]
                    and row['original'].get(f'chain-{side + 1}') == chain
                    for side, chain in enumerate(chains)
                )
                for row in original['fields']['contacts'][group]
            ),
            'mapping',
            'Both endpoints must retain their original participant correspondence.',
        )
        expected = distance_values(original['fields']['contacts'][group])
        quantity = QuantityRecord.from_dict(
            value['reported_distances'][group]
        ).to_quantity(
            field=distance_field(group),
            unit='angstrom',
            dimensionality={'[L]': 1},
            form='pint',
        )
        actual = np.asarray(puw.get_value(quantity, to_unit='angstrom'))
        _require(
            actual.shape == expected.shape
            and np.allclose(actual, expected, rtol=32 * np.finfo(float).eps, atol=0),
            'distances',
            'Preserve every selected original reported distance and its units.',
        )
    return payload(value)
