"""Native selection capture and data-only admission of its fixed membership.

MolSysMT alone evaluates selection syntax and molecular/spatial operations.
The inventory route supplies its resident topology/chemical domains without
structures. BindingMT checks the owned snapshot's scope and recorded axes.
"""

from numbers import Integral

import numpy as np

from bindingmt._private.interface_records import _require
from bindingmt._private.serialization import json_copy
from bindingmt._private.site_records import (
    capability,
    keys,
    membership_payload,
    scope_record,
    text,
)
from bindingmt._private.site_region import _index, _same

_COUNTS = ('n_atoms', 'n_groups', 'n_structures', 'n_chemical_states')
_REQUEST = {
    'selection',
    'element',
    'syntax',
    'structure_indices',
    'chemical_state',
    'max_source_atoms',
    'max_structure_occurrences',
    'max_dense_pairs',
}
_OBSERVATIONS = {
    *_COUNTS,
    'reference_chemical_state_index',
    'structure_chemical_state_index',
}


def _query(value, budget):
    if isinstance(value, str):
        _require(text(value), 'selection', 'Use a nonempty selection.')
        capability(
            '@' not in value, 'self-contained selection without external bindings'
        )
        return value
    if isinstance(value, Integral) and not isinstance(value, bool):
        return _index(value, 'selection')
    _require(
        isinstance(value, (list, tuple))
        or isinstance(value, np.ndarray)
        and value.ndim == 1,
        'selection',
        'Use a string or a scalar/flat integer selection.',
    )
    _require(
        len(value) <= budget, 'selection', 'Selection exceeds the admitted axis budget.'
    )
    values = [_index(i, 'selection') for i in value]
    _require(
        len(set(values)) == len(values),
        'selection',
        'A fixed membership needs unique entries; no implicit deduplication.',
    )
    return values


def _structures(value):
    if value is None:
        return None
    if isinstance(value, Integral) and not isinstance(value, bool):
        return [_index(value, 'structure_indices')]
    _require(
        isinstance(value, (list, tuple))
        or isinstance(value, np.ndarray)
        and value.ndim == 1,
        'structure_indices',
        'Use one explicit structure or None for inventory selection.',
    )
    _require(
        len(value) == 1,
        'structure_indices',
        'This profile selects exactly one structure; no aggregation.',
    )
    return [_index(value[0], 'structure_indices')]


def request_record(
    selection, element, syntax, structure_indices, chemical_state, **limits
):
    for field, limit in limits.items():
        _index(limit, field, minimum=1)
    capability(element in ('atom', 'group'), 'selection atom/group axis')
    capability(syntax == 'MolSysMT', 'native MolSysMT selection syntax')
    if chemical_state is not None and not (
        isinstance(chemical_state, str) and chemical_state in ('reference', 'structure')
    ):
        chemical_state = _index(chemical_state, 'chemical_state')
    return {
        'selection': _query(selection, limits['max_source_atoms']),
        'element': element,
        'syntax': syntax,
        'structure_indices': _structures(structure_indices),
        'chemical_state': chemical_state,
        **{k: int(v) for k, v in limits.items()},
    }


def _decision(request, observations):
    keys(request, _REQUEST, 'selection.request')
    admitted = request_record(**request)
    _require(
        _same(admitted, request), 'selection', 'Preserve canonical selection choices.'
    )
    keys(observations, _OBSERVATIONS, 'selection.observations')
    for name in _COUNTS:
        budget = (
            request['max_source_atoms']
            if name in ('n_atoms', 'n_groups')
            else request['max_structure_occurrences']
        )
        _index(observations[name], name, budget + 1)
    _require(
        observations['n_chemical_states'] > 0,
        'chemical_state',
        'Require a known native chemical-state collection.',
    )
    reference = observations['reference_chemical_state_index']
    associated = observations['structure_chemical_state_index']
    for state in (reference, associated):
        if state is not None:
            _index(state, 'chemical_state', observations['n_chemical_states'])
    structures = request['structure_indices']
    if structures is not None:
        _index(structures[0], 'structure_indices', observations['n_structures'])
        _require(
            associated is not None,
            'chemical_state',
            'Require the chosen structure association, without reference fallback.',
        )
    else:
        _require(
            associated is None,
            'chemical_state',
            'Inventory selection has no structure association.',
        )
    state = request['chemical_state']
    if state is None:
        state = 'reference' if structures is None else 'structure'
    if state == 'reference':
        selected, policy = reference, 'reference'
    elif state == 'structure':
        _require(
            structures is not None,
            'chemical_state',
            'Structure policy needs one explicit structure.',
        )
        selected, policy = associated, 'structure'
    else:
        selected, policy = state, 'explicit'
    selected = _index(selected, 'chemical_state', observations['n_chemical_states'])
    _require(
        structures is None or selected == associated,
        'chemical_state',
        'Explicit/reference state must agree with the chosen structure in this profile.',
    )
    # No local expression classifier. Any opaque structure-bound string could
    # invoke a full owner contact query and receives this conservative bound.
    pairs = (
        observations['n_atoms'] ** 2
        if structures is not None
        and isinstance(request['selection'], str)
        and request['selection'] != 'all'
        else 0
    )
    _require(
        pairs <= request['max_dense_pairs'],
        'max_dense_pairs',
        'Opaque structure-bound selection exceeds the full-source pair budget.',
    )
    return {
        'mode': 'inventory' if structures is None else 'single_structure',
        'chemical_state_policy': policy,
        'chemical_state_index': selected,
        'structure_indices': [] if structures is None else structures,
        'budgeted_dense_pairs': pairs,
    }


def validate_capture(payload, scope):
    """Check historical owned relationships without evaluating an expression."""
    capture = payload['selection']
    keys(
        capture,
        {
            'request',
            'effective',
            'observations',
            'provider',
            'membership_policy',
            'current_applicability',
        },
        'selection',
    )
    request, observations = capture['request'], capture['observations']
    effective = _decision(request, observations)
    _require(
        _same(effective, capture['effective']),
        'selection.effective',
        'Preserve effective captured context and resource choices.',
    )
    keys(capture['provider'], {'library', 'version'}, 'selection.provider')
    _require(
        capture['provider']['library'] == 'molsysmt'
        and text(capture['provider']['version'])
        and capture['membership_policy'] == 'fixed_snapshot'
        and capture['current_applicability'] == 'not_assessed',
        'selection',
        'Preserve historical owner provenance without asserting current applicability.',
    )
    expected = membership_payload(payload['source_indices'], payload['element'])
    expected.update(
        resolution={'state': 'resolved', 'cause': 'captured_selection'},
        selection=capture,
    )
    _require(
        _same(payload, expected) and payload['element'] == request['element'],
        'membership',
        'Preserve the evaluated membership and element axis.',
    )
    count = observations['n_atoms' if payload['element'] == 'atom' else 'n_groups']
    _require(
        len(payload['source_indices']) <= request['max_source_atoms'],
        'membership',
        'Result exceeds the admitted axis budget.',
    )
    for i in payload['source_indices']:
        _index(i, 'membership', count)
    query = request['selection']
    explicit = (
        (list(range(count)) if query == 'all' else None)
        if isinstance(query, str)
        else [query]
        if isinstance(query, int)
        else query
    )
    if explicit is not None:
        _require(
            _same(payload['source_indices'], explicit),
            'selection',
            'Preserve explicit numeric/all selection order without inferred members.',
        )
    structures = request['structure_indices']
    expected_scope = scope_record(
        scope['source_id'],
        structures,
        None if structures is None else [effective['chemical_state_index']],
    )
    _require(
        _same(scope, expected_scope),
        'scope',
        'Selection capture must retain its actual source structure/state scope.',
    )


def capture(molecular_system, request):
    import molsysmt as msm
    from molsysmt.native import MolSys

    capability(
        isinstance(molecular_system, MolSys),
        'resident native selection source',
    )

    def observe():
        counts = {
            name: _index(msm.get(molecular_system, **{name: True}), name)
            for name in _COUNTS
        }
        # Check axis budgets before state or selection-dependent operations.
        for name in _COUNTS:
            budget = (
                request['max_source_atoms']
                if name in ('n_atoms', 'n_groups')
                else request['max_structure_occurrences']
            )
            _index(counts[name], name, budget + 1)
        reference = msm.get(molecular_system, reference_chemical_state_index=True)
        association = None
        if request['structure_indices'] is not None:
            _index(
                request['structure_indices'][0],
                'structure_indices',
                counts['n_structures'],
            )
            state = msm.get(
                molecular_system,
                structure_indices=request['structure_indices'],
                structure_chemical_state_index=True,
            )[0]
            _require(
                isinstance(state, Integral) and not isinstance(state, bool),
                'chemical_state',
                'The chosen structure needs a known state association.',
            )
            association = int(state)
        return {
            **counts,
            'reference_chemical_state_index': None
            if reference is None
            else _index(reference, 'chemical_state'),
            'structure_chemical_state_index': association,
        }

    before = observe()
    effective = _decision(request, before)
    # Resident public owner domains; no extraction/copy/conversion and no local
    # parser. Without Structures the owner rejects coordinate-requiring queries.
    system = (
        [molecular_system.topology, molecular_system.chemical_states]
        if effective['mode'] == 'inventory'
        else molecular_system
    )

    def select():
        result = msm.select(
            system,
            selection=request['selection'],
            element=request['element'],
            syntax=request['syntax'],
            structure_indices=effective['structure_indices'],
            chemical_state=effective['chemical_state_index'],
        )
        _require(
            isinstance(result, (list, tuple, np.ndarray))
            and len(result) <= request['max_source_atoms'],
            'selection',
            'Require a bounded flat membership result.',
        )
        return membership_payload(result, request['element'])['source_indices']

    indices = select()
    _require(
        observe() == before,
        'molecular_system',
        'Observed selection context changed during capture.',
    )
    _require(
        select() == indices and observe() == before,
        'selection',
        'Selection changed during capture; supply consistent inputs.',
    )
    record = {
        'request': request,
        'effective': effective,
        'observations': before,
        'provider': {'library': 'molsysmt', 'version': msm.__version__},
        'membership_policy': 'fixed_snapshot',
        'current_applicability': 'not_assessed',
    }
    return indices, json_copy(record, 'selection')
