"""Native fixed-membership search boxes and bounded problem/use admission."""

import hashlib
import json
import re
from uuid import UUID, uuid4

import numpy as np
import pyunitwizard as puw

from bindingmt._private.interface_records import _digest, _require
from bindingmt._private.serialization import json_copy
from bindingmt._private.site_domain_geometry import _padding, _quantity, _sealed
from bindingmt._private.site_records import keys
from bindingmt._private.site_region import _index, _same
from bindingmt._private.smonitor import CapabilityMismatchError

FORMAT = 'bindingmt.site_domain_application'


def _supported(condition, name):
    if not condition:
        raise CapabilityMismatchError(capability=name, engine=FORMAT)


FIELDS = {
    name: f'{FORMAT}.{name}'
    for name in ('padding', 'minimum', 'maximum', 'center', 'lengths')
}
ASSURANCE = {
    'scope': 'placed_geometry',
    'placement': 'captured_source_positions',
    'current_finding': 'not_assessed',
    'chemical_equivalence': 'not_assessed',
    'source_association': 'caller_precondition',
    'lifetime': 'recheck_source_and_normalized_receptor_at_problem_construction_and_export',
    'consistent_inputs': 'caller_precondition',
}
LOSSES = {
    'support_relation': 'encloses_fixed_atom_positions',
    'added_space': 'unquantified',
    'binding_volume': 'not_inferred',
    'binding_support': 'unknown',
    'coordinates': 'cartesian_as_stored',
    'pbc': False,
}
BINDING = {
    'participant': 'receptor',
    'normalized_structure_index': 0,
    'normalization': 'native_source_sorted_unique_atom_subset_single_structure',
}


def _values(quantity):
    return np.asarray(puw.get_value(quantity, to_unit='nm'), dtype=float)


def _close(left, right):
    a, b = _values(left), _values(right)
    return (
        a.shape == b.shape == (3,)
        and np.isfinite(a).all()
        and np.isfinite(b).all()
        and np.allclose(a, b, rtol=8 * np.finfo(float).eps, atol=0)
    )


def _positions(system, atoms, structure):
    import molsysmt as msm

    quantity = msm.get(
        system,
        element='atom',
        selection=atoms,
        structure_indices=[structure],
        chemical_state='structure',
        coordinates=True,
    )
    values = np.array(puw.get_value(quantity, to_unit='nm'), dtype=float, copy=True)
    _require(
        values.shape == (1, len(atoms), 3) and np.isfinite(values).all(),
        'coordinates',
        'Require one finite Cartesian position per ordered support atom.',
    )
    values[values == 0] = 0.0  # A sign on zero is not positional information.
    return values[0]


def _position_digest(values):
    return hashlib.sha256(
        json.dumps(
            {'unit': 'nm', 'values': values.tolist()},
            sort_keys=True,
            separators=(',', ':'),
            allow_nan=False,
        ).encode()
    ).hexdigest()


def _box_from_bounds(record):
    from bindingmt.core.search_domain import BoxRegion

    minimum = _quantity(record['bounds']['minimum'], 'minimum', (3,), fields=FIELDS)
    maximum = _quantity(record['bounds']['maximum'], 'maximum', (3,), fields=FIELDS)
    padding = _quantity(
        record['request']['padding'], 'padding', (), sign='nonnegative', fields=FIELDS
    )
    return BoxRegion.from_bounds(
        puw.quantity(minimum, 'nm'),
        puw.quantity(maximum, 'nm'),
        padding=puw.quantity(float(padding), 'nm'),
    )


def validate_derivation(record, center, lengths):
    """Admit historical data; source authenticity/positions need actual use checks."""
    if (
        isinstance(record, dict)
        and record.get('format') == 'bindingmt.site_geometry_application'
    ):
        from bindingmt._private.site_geometry_application import (
            validate_derivation as validate,
        )

        return validate(record, center, lengths)
    from bindingmt.core.site import BindingSite

    record = json_copy(record, 'derivation')
    keys(
        record,
        {
            'format',
            'format_version',
            'experimental',
            'record_id',
            'region',
            'request',
            'positions',
            'bounds',
            'box',
            'assurance',
            'losses',
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
        'Preserve the explicit experimental application version and content digest.',
    )
    try:
        UUID(record['record_id'])
    except (ValueError, TypeError, AttributeError):
        _require(False, 'record_id', 'Use an application UUID.')
    keys(
        record['request'],
        {'method', 'padding', 'assurance_scope', 'context_policy'},
        'request',
    )
    _require(
        record['request']['method'] == 'enclosing_box'
        and record['request']['assurance_scope'] == 'placed_geometry'
        and record['request']['context_policy'] == 'require_match',
        'request',
        'Preserve the admitted enclosure and actual placement request.',
    )
    _require(
        _same(record['assurance'], ASSURANCE)
        and _same(record['losses'], LOSSES)
        and _same(record['consumer_binding'], BINDING),
        'assurance',
        'Preserve the actual guarantee, consumer boundary and missing stronger claims.',
    )
    region = BindingSite.region_from_dict(record['region'], experimental=True)
    _require(
        region['request']['derived'] is False
        and region['request']['coverage_policy'] == 'require_complete'
        and region['request']['context_policy'] == 'require_match'
        and region['projection']['coverage'] == 'complete'
        and bool(region['projection']['exported_target_atom_indices']),
        'region',
        'Require complete nonempty membership on the original native source axis.',
    )
    keys(record['positions'], {'encoding', 'sha256'}, 'positions')
    _require(
        record['positions']['encoding']
        == 'ordered_binary64_nm_json_signed_zero_normalized'
        and isinstance(record['positions']['sha256'], str)
        and re.fullmatch('[0-9a-f]{64}', record['positions']['sha256']) is not None,
        'positions',
        'Preserve the exact ordered-position witness; it authenticates no source.',
    )
    keys(record['bounds'], {'minimum', 'maximum'}, 'bounds')
    keys(record['box'], {'center', 'lengths'}, 'box')
    expected = _box_from_bounds(record)
    recorded_center = _quantity(record['box']['center'], 'center', (3,), fields=FIELDS)
    recorded_lengths = _quantity(
        record['box']['lengths'], 'lengths', (3,), sign='positive', fields=FIELDS
    )
    _require(
        np.array_equal(recorded_center, _values(expected.center))
        and np.array_equal(recorded_lengths, _values(expected.lengths)),
        'box',
        'Stored output must equal the retained bounds/padding recipe.',
    )
    _require(
        _close(center, expected.center) and _close(lengths, expected.lengths),
        'derivation',
        'Current BoxRegion geometry differs from its application recipe.',
    )
    keys(record['producer'], {'library', 'bindingmt_version', 'operation'}, 'producer')
    _require(
        record['producer']['library'] == 'bindingmt'
        and record['producer']['operation'] == 'to_search_domain'
        and isinstance(record['producer']['bindingmt_version'], str),
        'producer',
        'Preserve the original application producer.',
    )
    return record


def _current_region(record, source):
    from bindingmt.core.site import BindingSite

    original = record['region']
    site = BindingSite.from_dict(original['site'], experimental=True)
    request = original['request']
    current = site.resolve_region(
        source,
        **{
            k: request[k]
            for k in (
                'view',
                'structure_index',
                'target_structure_index',
                'chemical_state',
                'expansion',
                'coverage_policy',
                'context_policy',
                'max_source_atoms',
                'max_structure_occurrences',
            )
        },
    )
    _require(
        _same(current['projection'], original['projection'])
        and _same(current['observations'], original['observations']),
        'molecular_system',
        'Source membership/state/count witnesses changed since domain capture.',
    )
    return current


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
    assurance_scope,
    context_policy,
    max_source_atoms,
    max_structure_occurrences,
    frame_id=None,
):
    from bindingmt.core.search_domain import BoxRegion, _point_bounds
    from bindingmt.core.site import _producer

    if len(site.views) == 1 and site.views[0]['kind'] in ('sphere', 'box'):
        from bindingmt._private.site_geometry_application import derive as place

        return place(
            site,
            molecular_system,
            method=method,
            view=view,
            structure_index=structure_index,
            chemical_state=chemical_state,
            expansion=expansion,
            padding=padding,
            frame_id=frame_id,
            assurance_scope=assurance_scope,
            context_policy=context_policy,
            max_source_atoms=max_source_atoms,
            max_structure_occurrences=max_structure_occurrences,
        )
    _supported(
        frame_id is None,
        'membership enclosure without an independent declared geometry frame',
    )
    _supported(method == 'enclosing_box', 'native membership enclosure method')
    _supported(
        assurance_scope == 'placed_geometry',
        'placed geometry without current-finding assurance',
    )
    _supported(context_policy == 'require_match', 'checked native application context')
    _supported(
        len(site.views) == 1 and site.views[0]['kind'] == 'molecular_membership',
        'native atom/group view for molecular enclosure',
    )
    _require(
        padding is not None,
        'padding',
        'Declare explicit length padding for an enclosing projection.',
    )
    padding_record, padding_nm = _padding(padding, fields=FIELDS)
    region = site.resolve_region(
        molecular_system,
        view=view,
        structure_index=structure_index,
        chemical_state=chemical_state,
        expansion=expansion,
        coverage_policy='require_complete',
        context_policy=context_policy,
        max_source_atoms=max_source_atoms,
        max_structure_occurrences=max_structure_occurrences,
    )
    atoms = region['projection']['exported_target_atom_indices']
    _require(
        bool(atoms),
        'membership',
        'An empty fixed region defines no enclosing search box.',
    )
    frame = region['projection']['axes']['source_structure_index']
    values = _positions(molecular_system, atoms, frame)
    pad = puw.quantity(padding_nm, 'nm')
    minimum, maximum = _point_bounds(puw.quantity(values, 'nm'), pad)
    box = BoxRegion.from_bounds(
        minimum, maximum, padding=pad, name=site.to_dict()['name']
    )
    record = {
        'format': FORMAT,
        'format_version': '0.1',
        'experimental': True,
        'record_id': str(uuid4()),
        'region': region,
        'request': {
            'method': method,
            'padding': padding_record,
            'assurance_scope': assurance_scope,
            'context_policy': context_policy,
        },
        'positions': {
            'encoding': 'ordered_binary64_nm_json_signed_zero_normalized',
            'sha256': _position_digest(values),
        },
        'bounds': {
            'minimum': _sealed(_values(minimum), 'minimum', fields=FIELDS),
            'maximum': _sealed(_values(maximum), 'maximum', fields=FIELDS),
        },
        'box': {
            'center': _sealed(_values(box.center), 'center', fields=FIELDS),
            'lengths': _sealed(_values(box.lengths), 'lengths', fields=FIELDS),
        },
        'assurance': dict(ASSURANCE),
        'losses': dict(LOSSES),
        'consumer_binding': dict(BINDING),
        'producer': _producer('to_search_domain'),
    }
    _current_region(record, molecular_system)
    _require(
        _position_digest(_positions(molecular_system, atoms, frame))
        == record['positions']['sha256'],
        'molecular_system',
        'Source positions changed during enclosure capture.',
    )
    record = json_copy(record, 'derivation')
    record['content_sha256'] = _digest(record)
    box._set_derivation(record)
    return box


def validate_problem_application(problem, domain, *, checkpoint='problem_construction'):
    """Recheck actual source and normalized receptor immediately before use."""
    recipe = domain.derivation
    if (
        isinstance(recipe, dict)
        and recipe.get('format') == 'bindingmt.site_geometry_application'
    ):
        from bindingmt._private.site_geometry_application import (
            validate_problem_application as validate,
        )

        return validate(problem, domain, checkpoint=checkpoint)
    import molsysmt as msm
    from molsysmt.native import MolSys

    from bindingmt.core.problem import DockingProblem
    from bindingmt.core.search_domain import BoxRegion

    _require(
        isinstance(problem, DockingProblem),
        'problem',
        'A placed site box needs its actual DockingProblem at export.',
    )
    record = domain.derivation
    _require(
        isinstance(record, dict) and record.get('format') == FORMAT,
        'derivation',
        'Require a native site application.',
    )
    construction = problem.search_domain_application
    if construction is not None:
        _require(
            construction['application_ref']
            == {
                'record_id': record['record_id'],
                'content_sha256': record['content_sha256'],
            },
            'search_domain_application',
            'A replacement application requires a newly checked problem.',
        )
    _require(
        isinstance(problem.search_domain, BoxRegion)
        and _same(problem.search_domain.derivation, record),
        'problem',
        'The exported domain must be the application assigned to this problem.',
    )
    _supported(
        isinstance(problem.receptor, MolSys),
        'original resident native receptor',
    )
    current = _current_region(record, problem.receptor)
    projection = current['projection']
    atoms = projection['exported_target_atom_indices']
    axes = projection['axes']
    _require(
        problem.receptor_structure_index == axes['source_structure_index'],
        'problem',
        'Problem receptor must choose the original application structure.',
    )
    source_atoms = problem.receptor_atom_indices
    count = current['observations']['source']['n_atoms']
    _require(
        isinstance(source_atoms, list)
        and bool(source_atoms)
        and source_atoms == sorted(set(source_atoms)),
        'problem',
        'Require the admitted sorted unique native extraction axis.',
    )
    for atom in source_atoms:
        _index(atom, 'problem.receptor_atom_indices', count)
    source_to_local = {atom: local for local, atom in enumerate(source_atoms)}
    _require(
        all(atom in source_to_local for atom in atoms),
        'problem',
        'Every fixed support atom must survive receptor selection.',
    )
    local_atoms = [source_to_local[atom] for atom in atoms]
    normalized = problem.receptor_molsys
    _supported(
        isinstance(normalized, MolSys),
        'actual native normalized receptor',
    )
    _require(
        msm.get(normalized, n_structures=True) == 1
        and msm.get(normalized, n_atoms=True) == len(source_atoms),
        'problem',
        'Normalized receptor axis counts contradict its source relation.',
    )
    state = msm.get(
        normalized, structure_indices=[0], structure_chemical_state_index=True
    )[0]
    state = _index(
        state, 'chemical_state', current['observations']['source']['n_chemical_states']
    )
    _require(
        state == axes['chemical_state_index']
        and problem.receptor_chemical_state_index == state,
        'chemical_state',
        'Actual normalized association must preserve the selected source state position.',
    )
    original_positions = _positions(
        problem.receptor, atoms, axes['source_structure_index']
    )
    normalized_positions = _positions(normalized, local_atoms, 0)
    digest = record['positions']['sha256']
    _require(
        _position_digest(original_positions) == digest
        and _position_digest(normalized_positions) == digest,
        'coordinates',
        'Original or actual normalized support positions changed; recapture the application.',
    )
    padding = _quantity(
        record['request']['padding'], 'padding', (), sign='nonnegative', fields=FIELDS
    )
    expected = BoxRegion.from_points(
        puw.quantity(normalized_positions, 'nm'),
        padding=puw.quantity(float(padding), 'nm'),
    )
    _require(
        _close(domain.center, expected.center)
        and _close(domain.lengths, expected.lengths),
        'search_domain',
        'Application box does not enclose the actual normalized input by its recorded rule.',
    )
    return {
        'format': 'bindingmt.site_domain_use',
        'format_version': '0.1',
        'experimental': True,
        'application_ref': {
            'record_id': record['record_id'],
            'content_sha256': record['content_sha256'],
        },
        'checked_at': checkpoint,
        'consumer': 'DockingProblem.receptor',
        'source_structure_index': axes['source_structure_index'],
        'normalized_structure_index': 0,
        'chemical_state_index': state,
        'normalized_atom_source_indices': source_atoms,
        'support_source_atom_indices': atoms,
        'support_local_atom_indices': local_atoms,
        'positions_sha256': digest,
        'lifetime': 'historical construction receipt; export requires a new actual-input check',
    }


def validate_staged_receptor(problem, prepared, exported, payload, backend_box):
    """Check owner-mapped prepared positions and the exact rounded engine input.

    The owner writer/report binds the supplied bytes; its public file handler
    reads positions. This is a Cartesian box admission, not chemical validation.
    """
    import molsysmt as msm
    from molsysmt.native import PDBQTFileHandler

    from bindingmt.core.search_domain import _descriptor
    from bindingmt.preparation import PreparedReceptor

    use = validate_problem_application(
        problem, problem.search_domain, checkpoint='backend_export'
    )
    recipe = problem.search_domain.derivation
    geometric = recipe['format'] == 'bindingmt.site_geometry_application'
    _require(
        isinstance(prepared, PreparedReceptor),
        'receptor',
        'Require the actual native prepared receptor export.',
    )
    limit = (
        recipe['request']['max_source_atoms']
        if geometric
        else recipe['region']['request']['max_source_atoms']
    )
    _require(
        0 < prepared.n_atoms <= limit,
        'receptor',
        'Prepared axis exceeds the admitted native application budget.',
    )
    projection = prepared.metadata.get('molecular_projection')
    _require(
        isinstance(projection, dict),
        'receptor',
        'Native preparation must retain its owner atom-map report.',
    )
    _require(
        payload == exported.text.encode('utf-8'),
        'receptor',
        'Staged bytes differ from the actual owner export.',
    )
    payload_sha = hashlib.sha256(payload).hexdigest()
    _require(
        exported.atom_map.get('payload_sha256') == payload_sha,
        'receptor',
        'Owner export map must bind the exact submitted bytes.',
    )
    # These are actual bounded operator reports, not inferred label matches.
    normalized_count = msm.get(problem.receptor_molsys, n_atoms=True)
    _require(
        projection.get('n_source_atoms') == normalized_count
        and projection.get('n_output_atoms') == prepared.n_atoms
        and exported.atom_map.get('n_source_atoms') == prepared.n_atoms
        and exported.atom_map.get('n_output_atoms') == prepared.n_atoms
        and 'stages' not in projection
        and 'stages' not in exported.atom_map,
        'receptor',
        'Require bounded atomic preparation/export reports before map normalization.',
    )
    preparation_map = msm.compose_atom_maps([projection])
    _require(
        preparation_map['n_source_atoms'] == normalized_count
        and preparation_map['n_output_atoms'] == prepared.n_atoms,
        'receptor',
        'Preparation report does not cover the normalized receptor axes.',
    )
    combined = msm.compose_atom_maps([projection, exported.atom_map])
    _require(
        combined['n_output_atoms'] == prepared.n_atoms,
        'receptor',
        'Writer axis must match the prepared atom count.',
    )
    checked_local = (
        use['checked_local_atom_indices']
        if geometric
        else use['support_local_atom_indices']
    )
    checked_source = (
        use['checked_source_atom_indices']
        if geometric
        else use['support_source_atom_indices']
    )
    prepared_atoms = [
        int(preparation_map['source_to_output_atom_indices'][i]) for i in checked_local
    ]
    written_atoms = [
        int(combined['source_to_output_atom_indices'][i]) for i in checked_local
    ]
    _require(
        all(i >= 0 for i in prepared_atoms + written_atoms),
        'receptor',
        'Preparation omitted a checked original atom; direct frame identity is required.'
        if geometric
        else 'Preparation omitted a fixed support atom; representative identity is insufficient.',
    )
    native_positions = _positions(exported.molecular_system, prepared_atoms, 0)
    _require(
        _position_digest(native_positions) == use['positions_sha256'],
        'receptor',
        'Actual prepared full-precision support positions differ from the checked application.',
    )
    with PDBQTFileHandler('pdbqt_text:' + payload.decode('utf-8')) as reader:
        _require(
            reader.n_structures == 1 and reader.n_atoms == prepared.n_atoms,
            'receptor',
            'Actual PDBQT axes contradict the export report.',
        )
        points = np.asarray(
            puw.get_value(
                reader.get_coordinates(indices=written_atoms, structure_indices=[0]),
                to_unit='angstrom',
            )
        )
    _require(
        points.shape == (1, len(written_atoms), 3) and np.isfinite(points).all(),
        'receptor',
        'Require finite actual written support positions.',
    )
    center, size = np.asarray(backend_box['center']), np.asarray(backend_box['size'])
    _require(
        backend_box['unit'] == 'angstrom'
        and center.shape == size.shape == (3,)
        and np.isfinite(center).all()
        and np.isfinite(size).all()
        and (size > 0).all(),
        'backend_box',
        'Require the actual finite positive numeric engine box.',
    )
    expected_center = [
        round(float(v), 6)
        for v in puw.get_value(problem.search_domain.center, to_unit='angstrom')
    ]
    expected_size = [
        round(float(v), 6)
        for v in puw.get_value(problem.search_domain.lengths, to_unit='angstrom')
    ]
    _require(
        _same(backend_box['center'], expected_center)
        and _same(backend_box['size'], expected_size),
        'backend_box',
        'The actual rounded engine box differs from the checked application.',
    )
    if geometric:
        from bindingmt._private.site_geometry_application import validate_numeric_box

        expected_written = np.asarray(
            [
                [float(format(float(v), '.3f')) for v in row]
                for row in puw.get_value(
                    puw.quantity(native_positions, 'nm'), to_unit='angstrom'
                )
            ]
        )
        _require(
            np.allclose(
                points[0], expected_written, rtol=8 * np.finfo(float).eps, atol=0
            ),
            'receptor',
            'Written receptor positions differ from the checked frame at owner PDBQT precision.',
        )
        validate_numeric_box(recipe, center, size)
        return {
            'format': 'bindingmt.site_geometry_payload_use',
            'format_version': '0.1',
            'experimental': True,
            'application_ref': use['application_ref'],
            'receptor_payload_sha256': payload_sha,
            'normalized_to_pdbqt_atom_map': _descriptor(combined),
            'checked_source_atom_indices': checked_source,
            'checked_normalized_atom_indices': checked_local,
            'checked_pdbqt_atom_indices': written_atoms,
            'backend_box': json_copy(backend_box, 'backend_box'),
            'guarantee': 'checked original frame positions; effective box encloses declared geometry',
            'frame_association': 'caller_precondition',
            'written_coordinate_policy': 'owner three-decimal angstrom formatting with binary64 representation roundoff',
            'membership': 'not_requested',
            'chemical_equivalence': 'not_assessed',
            'current_finding': 'not_assessed',
        }
    _require(
        (points[0] >= center - size / 2).all()
        and (points[0] <= center + size / 2).all(),
        'backend_box',
        'Written PDBQT support is outside the actual rounded engine box; declare adequate padding.',
    )
    return {
        'format': 'bindingmt.site_domain_payload_use',
        'format_version': '0.1',
        'experimental': True,
        'application_ref': use['application_ref'],
        'receptor_payload_sha256': payload_sha,
        'normalized_to_pdbqt_atom_map': _descriptor(combined),
        'support_source_atom_indices': use['support_source_atom_indices'],
        'support_normalized_atom_indices': use['support_local_atom_indices'],
        'support_pdbqt_atom_indices': written_atoms,
        'backend_box': json_copy(backend_box, 'backend_box'),
        'guarantee': 'retained full-precision source positions; written support inside actual rounded engine box',
        'chemical_equivalence': 'not_assessed',
        'current_finding': 'not_assessed',
    }
