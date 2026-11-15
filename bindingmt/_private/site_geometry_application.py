"""Caller-associated declared geometry with checked native frame/use context."""

import re
from uuid import UUID, uuid4

import pyunitwizard as puw

from bindingmt._private.interface_records import _digest, _require
from bindingmt._private.serialization import json_copy
from bindingmt._private.site_application import BINDING, _position_digest, _positions
from bindingmt._private.site_domain_geometry import (
    _geometry,
    derive_geometry_box,
)
from bindingmt._private.site_domain_geometry import (
    validate_derivation as validate_offline,
)
from bindingmt._private.site_records import keys, text
from bindingmt._private.site_region import _index, _same
from bindingmt._private.smonitor import CapabilityMismatchError

FORMAT = 'bindingmt.site_geometry_application'
ASSURANCE = {
    'scope': 'placed_geometry',
    'placement': 'caller_associated_source_frame_with_checked_native_positions',
    'frame_association': 'caller_precondition',
    'source_association': 'caller_precondition',
    'consistent_inputs': 'caller_precondition',
    'current_finding': 'not_assessed',
    'chemical_equivalence': 'not_assessed',
    'membership': 'not_requested',
    'lifetime': 'recheck_whole_source_frame_and_normalized_receptor_at_construction_and_export',
}
_COUNTS = ('n_atoms', 'n_groups', 'n_structures', 'n_chemical_states')
_REQUEST = {
    'frame_id',
    'structure_index',
    'chemical_state',
    'assurance_scope',
    'context_policy',
    'max_source_atoms',
    'max_structure_occurrences',
}


def _supported(condition, capability):
    if not condition:
        raise CapabilityMismatchError(capability=capability, engine=FORMAT)


def _choose(site, request, counts):
    for field in ('max_source_atoms', 'max_structure_occurrences'):
        _index(request[field], field, minimum=1)
    for field in _COUNTS:
        limit = (
            request['max_source_atoms']
            if field in ('n_atoms', 'n_groups')
            else request['max_structure_occurrences']
        )
        _index(counts[field], field, limit + 1)
    _require(
        all(counts[k] > 0 for k in ('n_atoms', 'n_structures', 'n_chemical_states')),
        'molecular_system',
        'Require a nonempty native frame and known chemical state.',
    )
    view = site.views[0]
    scope = view['scope']
    declared = scope.get('source_structure_indices', [])
    for i in declared:
        _index(i, 'scope', counts['n_structures'])
    chosen = request['structure_index']
    if chosen is None:
        _require(
            len(declared) == 1 or not declared and counts['n_structures'] == 1,
            'structure_index',
            'Choose one original structure explicitly.',
        )
        chosen = declared[0] if declared else 0
    chosen = _index(chosen, 'structure_index', counts['n_structures'])
    _require(
        not declared or chosen in declared,
        'structure_index',
        'The source structure contradicts the declared geometry scope.',
    )
    state = request['chemical_state']
    if state is not None and not (isinstance(state, str) and state == 'structure'):
        _index(state, 'chemical_state', counts['n_chemical_states'])
    return chosen


def _check_context(site, request, observation):
    chosen = _choose(site, request, observation)
    actual = _index(
        observation['structure_index'], 'structure_index', observation['n_structures']
    )
    _require(
        actual == chosen,
        'structure_index',
        'Preserve the actual selected original frame.',
    )
    state = _index(
        observation['chemical_state_index'],
        'chemical_state',
        observation['n_chemical_states'],
    )
    requested = request['chemical_state']
    if requested is not None and not (
        isinstance(requested, str) and requested == 'structure'
    ):
        _require(
            state == requested,
            'chemical_state',
            'The requested state contradicts the structure association.',
        )
    scope = site.views[0]['scope']
    declared_states = scope.get('chemical_state_indices', [])
    if declared_states and declared_states[0] is not None:
        _require(
            declared_states == [state],
            'chemical_state',
            "Preserve the geometry declaration's known state position.",
        )


def _observe(site, request, system):
    import molsysmt as msm
    from molsysmt.native import MolSys

    _supported(
        isinstance(system, MolSys),
        'original resident native geometry source',
    )
    counts = {k: _index(msm.get(system, **{k: True}), k) for k in _COUNTS}
    structure = _choose(site, request, counts)
    state = msm.get(
        system, structure_indices=[structure], structure_chemical_state_index=True
    )[0]
    result = {
        **counts,
        'structure_index': structure,
        'chemical_state_index': _index(
            state, 'chemical_state', counts['n_chemical_states']
        ),
    }
    _check_context(site, request, result)
    return result


def _source_site(record):
    from bindingmt.core.site import BindingSite

    return BindingSite.from_dict(record['geometry']['site'], experimental=True)


def _current(record, source):
    site = _source_site(record)
    before = _observe(site, record['request'], source)
    _require(
        _same(before, record['observation']),
        'molecular_system',
        'Original source counts or selected frame/state changed; recapture placement.',
    )
    points = _positions(
        source, list(range(before['n_atoms'])), before['structure_index']
    )
    _require(
        _position_digest(points) == record['positions_sha256'],
        'coordinates',
        'Captured source frame changed; explicitly recapture geometry placement.',
    )
    _require(
        _same(before, _observe(site, record['request'], source)),
        'molecular_system',
        'Source context changed during checking; supply consistent inputs.',
    )
    return points


def derive(
    site,
    molecular_system,
    *,
    method,
    view,
    structure_index,
    chemical_state,
    expansion,
    padding,
    frame_id,
    assurance_scope,
    context_policy,
    max_source_atoms,
    max_structure_occurrences,
):
    from bindingmt.core.site import _producer

    _supported(
        assurance_scope == 'placed_geometry',
        'placed geometry without current-finding assurance',
    )
    _supported(
        context_policy == 'require_match', 'checked native geometry placement context'
    )
    _supported(
        expansion is None, 'geometry placement without inferred molecular membership'
    )
    box = derive_geometry_box(site, method=method, view=view, padding=padding)
    declared_frame = site.views[0]['payload']['frame']
    _require(
        text(frame_id) and frame_id == declared_frame['id'],
        'frame_id',
        "Explicitly associate the saved frame ID with this source's stored Cartesian frame.",
    )
    request = {
        'frame_id': frame_id,
        'structure_index': structure_index,
        'chemical_state': chemical_state,
        'assurance_scope': assurance_scope,
        'context_policy': context_policy,
        'max_source_atoms': max_source_atoms,
        'max_structure_occurrences': max_structure_occurrences,
    }
    observation = _observe(site, request, molecular_system)
    points = _positions(
        molecular_system,
        list(range(observation['n_atoms'])),
        observation['structure_index'],
    )
    record = {
        'format': FORMAT,
        'format_version': '0.1',
        'experimental': True,
        'record_id': str(uuid4()),
        'geometry': box.derivation,
        'request': request,
        'observation': observation,
        'positions_sha256': _position_digest(points),
        'assurance': dict(ASSURANCE),
        'consumer_binding': dict(BINDING),
        'producer': _producer('to_search_domain'),
    }
    _current(record, molecular_system)
    record = json_copy(record, 'derivation')
    record['content_sha256'] = _digest(record)
    box._set_derivation(record)
    return box


def validate_derivation(record, center, lengths):
    record = json_copy(record, 'derivation')
    keys(
        record,
        {
            'format',
            'format_version',
            'experimental',
            'record_id',
            'geometry',
            'request',
            'observation',
            'positions_sha256',
            'assurance',
            'consumer_binding',
            'producer',
            'content_sha256',
        },
        'derivation',
    )
    _require(
        record['format'] == FORMAT
        and record['format_version'] == '0.1'
        and record['experimental'] is True
        and record['content_sha256'] == _digest(record),
        'derivation',
        'Preserve the experimental geometry application and content digest.',
    )
    try:
        UUID(record['record_id'])
    except (TypeError, ValueError, AttributeError):
        _require(False, 'record_id', 'Use a geometry application UUID.')
    validate_offline(record['geometry'], center, lengths)
    site = _source_site(record)
    keys(record['request'], _REQUEST, 'request')
    request = record['request']
    _require(
        text(request['frame_id'])
        and request['frame_id'] == site.views[0]['payload']['frame']['id']
        and request['assurance_scope'] == 'placed_geometry'
        and request['context_policy'] == 'require_match',
        'request',
        'Preserve explicit caller frame association and actual checked-context request.',
    )
    keys(
        record['observation'],
        {*_COUNTS, 'structure_index', 'chemical_state_index'},
        'observation',
    )
    _check_context(site, request, record['observation'])
    _require(
        isinstance(record['positions_sha256'], str)
        and re.fullmatch('[0-9a-f]{64}', record['positions_sha256']) is not None,
        'positions_sha256',
        'Preserve the exact whole-frame numerical witness.',
    )
    _require(
        _same(record['assurance'], ASSURANCE)
        and _same(record['consumer_binding'], BINDING),
        'assurance',
        'Preserve checked facts, caller preconditions and unsupported stronger claims.',
    )
    keys(record['producer'], {'library', 'bindingmt_version', 'operation'}, 'producer')
    _require(
        record['producer']['library'] == 'bindingmt'
        and isinstance(record['producer']['bindingmt_version'], str)
        and record['producer']['operation'] == 'to_search_domain',
        'producer',
        'Preserve the original producer.',
    )
    return record


def validate_problem_application(problem, domain, *, checkpoint):
    import molsysmt as msm
    from molsysmt.native import MolSys

    from bindingmt.core.problem import DockingProblem
    from bindingmt.core.search_domain import BoxRegion

    _require(
        isinstance(problem, DockingProblem),
        'problem',
        'Placed geometry export needs its actual DockingProblem.',
    )
    record = domain.derivation
    _require(
        isinstance(problem.search_domain, BoxRegion)
        and _same(problem.search_domain.derivation, record),
        'problem',
        'Export the actual application assigned to this problem.',
    )
    ref = {'record_id': record['record_id'], 'content_sha256': record['content_sha256']}
    previous = problem.search_domain_application
    _require(
        previous is None or _same(previous['application_ref'], ref),
        'search_domain_application',
        'A replacement application requires a newly checked problem.',
    )
    original = _current(record, problem.receptor)
    observation = record['observation']
    frame = observation['structure_index']
    _require(
        problem.receptor_structure_index == frame,
        'problem',
        'Problem must choose the captured source structure.',
    )
    atoms = problem.receptor_atom_indices
    _require(
        isinstance(atoms, list) and bool(atoms) and atoms == sorted(set(atoms)),
        'problem',
        'Require the known sorted unique native extraction axis.',
    )
    for atom in atoms:
        _index(atom, 'problem.receptor_atom_indices', observation['n_atoms'])
    normalized = problem.receptor_molsys
    _supported(
        isinstance(normalized, MolSys),
        'actual normalized native receptor',
    )
    _require(
        msm.get(normalized, n_structures=True) == 1
        and msm.get(normalized, n_atoms=True) == len(atoms),
        'problem',
        'Actual normalized axes disagree with their source relation.',
    )
    state = msm.get(
        normalized, structure_indices=[0], structure_chemical_state_index=True
    )[0]
    state = _index(state, 'chemical_state', observation['n_chemical_states'])
    _require(
        state == observation['chemical_state_index']
        and state == problem.receptor_chemical_state_index,
        'chemical_state',
        'Actual normalized state contradicts the original association.',
    )
    local = list(range(len(atoms)))
    points = _positions(normalized, local, 0)
    digest = _position_digest(original[atoms])
    _require(
        _position_digest(points) == digest,
        'coordinates',
        'Actual normalized receptor positions differ from the captured original frame.',
    )
    return {
        'format': 'bindingmt.site_geometry_use',
        'format_version': '0.1',
        'experimental': True,
        'application_ref': ref,
        'checked_at': checkpoint,
        'consumer': 'DockingProblem.receptor',
        'source_structure_index': frame,
        'normalized_structure_index': 0,
        'chemical_state_index': state,
        'normalized_atom_source_indices': atoms,
        'checked_source_atom_indices': atoms,
        'checked_local_atom_indices': local,
        'positions_sha256': digest,
        'source_frame_positions_sha256': record['positions_sha256'],
        'frame_id': record['request']['frame_id'],
        'lifetime': 'historical construction receipt; actual export rechecks captured frame and normalized input',
    }


def validate_numeric_box(record, center, size):
    """The effective rounded box must still enclose original declared support."""
    site = _source_site(record)
    request = record['geometry']['request']
    unpadded_center, unpadded_lengths, _ = _geometry(
        site, request['method'], request['view'], 0.0
    )
    lower = puw.get_value(
        puw.quantity(unpadded_center - unpadded_lengths / 2, 'nm'), to_unit='angstrom'
    )
    upper = puw.get_value(
        puw.quantity(unpadded_center + unpadded_lengths / 2, 'nm'), to_unit='angstrom'
    )
    _require(
        (center - size / 2 <= lower).all() and (center + size / 2 >= upper).all(),
        'backend_box',
        'Rounded engine box loses declared geometry; provide adequate explicit padding.',
    )
