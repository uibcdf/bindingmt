"""Actual native application checks original source and normalized engine input."""

import errno
import json

import molsysmt as msm
import numpy as np
import pandas as pd
import pytest
import pyunitwizard as puw
from molsysmt.native import MolSys

import bindingmt as bmt
from bindingmt._private.interface_records import _digest
from bindingmt._private.site_application import FIELDS
from bindingmt._private.site_domain_geometry import _sealed
from bindingmt._private.smonitor import ArgumentError, CapabilityMismatchError
from bmt_test_support.pdbqt_fixtures import MINIMAL_LIG_PDBQT

SOURCE = {'owner': 'authored.control', 'id': 'application-source', 'revision': '1'}


@pytest.fixture
def source():
    system = MolSys(n_atoms=6, n_groups=2, n_chains=1, n_molecules=2, n_entities=2)
    top = system.topology
    top.atoms['atom_id'] = ['same'] * 6
    top.atoms['atom_name'] = list('abcdef')
    top.atoms['atom_type'] = ['C'] * 6
    top.atoms['group_index'] = pd.array([0, 0, 0, 0, 0, 1], dtype='Int64')
    top.atoms['chain_index'] = pd.array([0] * 6, dtype='Int64')
    top.groups['group_id'] = ['same', 'same']
    top.groups['group_name'] = ['LIG', 'LIG']
    top.groups['group_type'] = ['small molecule'] * 2
    top.groups['molecule_index'] = pd.array([0, 1], dtype='Int64')
    top.molecules['entity_index'] = pd.array([0, 1], dtype='Int64')
    top.chains['chain_id'] = ['A']
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
            [0.25, 0.5, 1.0],
            [8.0, 8.0, 8.0],
        ]
    )
    system.structures.coordinates = puw.quantity(
        np.stack([points, points + [10.0, 0.0, 0.0], points + [4.0, 1.0, -2.0]]), 'nm'
    )
    system.structures.time = puw.quantity([0.0, 1.0, 2.0], 'ps')
    system.chemical_states.append_state()
    msm.set(system, chemical_state=0, formal_charge=[1, 0, 0, 0, 0, 0])
    msm.set(system, chemical_state=1, formal_charge=[0, 1, 0, 0, 0, 0])
    msm.set(system, structure_chemical_state_index=[0, 1, 1])
    return system


def declare(indices=(4, 1, 3)):
    return bmt.BindingSite.from_atoms(
        source=SOURCE,
        atom_indices=list(indices),
        structure_indices=[2],
        chemical_state_indices=[1],
    )


def domain(source, indices=(4, 1, 3), **options):
    return declare(indices).to_search_domain(
        source, method='enclosing_box', padding='0.2 nm', **options
    )


def problem(source, box, **options):
    return bmt.DockingProblem(
        source, MINIMAL_LIG_PDBQT, box, receptor_structure_index=2, **options
    )


def test_source2_to_local0_and_nonidentity_support_map(source):
    from bindingmt.engines.vina import _vina_box

    assert declare().inspect()['capabilities']['domain_placement'] is True
    assert declare().inspect()['molecular_resolution'] == 'not_performed'
    box = domain(source)
    p = problem(source, box, receptor_selection=[4, 3, 1])
    receipt = p.search_domain_application
    assert (
        receipt['source_structure_index'] == 2
        and receipt['normalized_structure_index'] == 0
    )
    assert receipt['normalized_atom_source_indices'] == [1, 3, 4]
    assert receipt['support_source_atom_indices'] == [4, 1, 3]
    assert receipt['support_local_atom_indices'] == [2, 0, 1]
    assert receipt['chemical_state_index'] == 1
    np.testing.assert_allclose(
        puw.get_value(box.center, to_unit='nm'), [4.5, 1.25, -0.5]
    )
    np.testing.assert_allclose(
        puw.get_value(box.lengths, to_unit='nm'), [1.4, 0.9, 3.4]
    )
    exported = box.to_backend_box(problem=p)
    np.testing.assert_allclose(exported['center'], [45, 12.5, -5])
    np.testing.assert_allclose(exported['lengths'], [14, 9, 34])
    assert _vina_box(box, problem=p) == ([45.0, 12.5, -5.0], [14.0, 9.0, 34.0])
    with pytest.raises(ArgumentError):
        box.to_backend_box()


@pytest.mark.parametrize('form', ['pint', 'unyt'])
def test_application_active_form_and_unit_policy(source, form):
    with puw.context(default_form=form, standard_units=['pm', 'fs']):
        before = puw.configure.report()
        box = domain(source)
        p = problem(source, box)
        assert (
            bmt.BoxRegion.from_dict(
                json.loads(json.dumps(box.to_dict(), sort_keys=True))
            ).to_dict()
            == box.to_dict()
        )
        assert box.to_backend_box(unit='nm', problem=p)['unit'] == 'nm'
        assert puw.configure.report() == before


def test_problem_roundtrip_keeps_application_and_rechecks_actual_inputs(
    source, tmp_path
):
    box = domain(source)
    p = problem(source, box)
    record = p.to_dict()
    assert record['search_domain_application'] == p.search_domain_application
    restored = bmt.DockingProblem.from_dict(
        record, receptor=source, partner=MINIMAL_LIG_PDBQT
    )
    assert restored.to_dict() == record
    (tmp_path / 'application-box.json').write_text(
        json.dumps(box.to_dict(), sort_keys=True)
    )
    (tmp_path / 'application-problem.json').write_text(
        json.dumps(record, sort_keys=True)
    )


def test_support_changes_refuse_on_source_and_normalized_input(source):
    box = domain(source)
    p = problem(source, box)
    msm.set(
        p.receptor_molsys,
        selection=[1],
        structure_indices=[0],
        coordinates=puw.quantity([[[99.0, 0.0, 0.0]]], 'nm'),
    )
    with pytest.raises(ArgumentError):
        box.to_backend_box(problem=p)
    # Historical geometry remains inspectable; export is the live-use boundary.
    assert bmt.BoxRegion.from_dict(box.to_dict()).to_dict() == box.to_dict()


@pytest.mark.parametrize('location', ['source', 'normalized'])
def test_interior_position_change_refuses_even_when_box_bounds_are_equal(
    source, location
):
    box = domain(source, indices=[4, 0, 3, 1, 2])
    p = problem(source, box)
    system = source if location == 'source' else p.receptor_molsys
    frame = 2 if location == 'source' else 0
    msm.set(
        system,
        selection=[4],
        structure_indices=[frame],
        coordinates=puw.quantity([[[4.4, 1.6, -0.9]]], 'nm'),
    )
    coords = msm.get(
        system, selection=[4, 0, 3, 1, 2], structure_indices=[frame], coordinates=True
    )
    equivalent_bounds = bmt.BoxRegion.from_points(coords[0], padding='0.2 nm')
    np.testing.assert_allclose(
        puw.get_value(equivalent_bounds.center, to_unit='nm'),
        puw.get_value(box.center, to_unit='nm'),
    )
    np.testing.assert_allclose(
        puw.get_value(equivalent_bounds.lengths, to_unit='nm'),
        puw.get_value(box.lengths, to_unit='nm'),
    )
    with pytest.raises(ArgumentError):
        box.to_backend_box(problem=p)
    assert p.to_dict()['search_domain_application'] == p.search_domain_application


def test_change_outside_fixed_support_does_not_claim_current_finding(source):
    box = domain(source)
    p = problem(source, box)
    msm.set(
        source,
        selection=[5],
        structure_indices=[2],
        coordinates=puw.quantity([[[100.0, 0.0, 0.0]]], 'nm'),
    )
    msm.set(
        p.receptor_molsys,
        selection=[5],
        structure_indices=[0],
        coordinates=puw.quantity([[[200.0, 0.0, 0.0]]], 'nm'),
    )
    assert box.to_backend_box(problem=p)['unit'] == 'angstrom'
    assert box.derivation['assurance']['current_finding'] == 'not_assessed'
    assert bmt.BoxRegion.from_dict(box.to_dict()).to_dict() == box.to_dict()


@pytest.mark.parametrize(
    'options',
    [
        {'receptor_selection': [1, 3]},
        {'receptor_structure_index': 1},
    ],
)
def test_problem_requires_all_members_and_original_structure(source, options):
    params = {'receptor_structure_index': 2, **options}
    with pytest.raises(ArgumentError):
        bmt.DockingProblem(source, MINIMAL_LIG_PDBQT, domain(source), **params)


@pytest.mark.parametrize('element', ['atom', 'group', 'selected'])
def test_supported_membership_profiles_propagate(source, element):
    if element == 'group':
        site = bmt.BindingSite.from_groups(
            source=SOURCE,
            group_indices=[0],
            structure_indices=[2],
            chemical_state_indices=[1],
        )
    elif element == 'selected':
        site = bmt.BindingSite.from_selection(
            source, [4, 1, 3], source=SOURCE, structure_indices=[2]
        )
    else:
        site = declare()
    box = site.to_search_domain(
        source,
        method='enclosing_box',
        padding='2 angstrom',
        expansion='source_members' if element == 'group' else None,
    )
    p = problem(source, box)
    assert box.derivation['region']['site'] == site.to_dict()
    assert box.to_backend_box(problem=p)['unit'] == 'angstrom'


@pytest.mark.parametrize(
    'options',
    [
        {'method': 'sphere_enclosure'},
        {'assurance_scope': 'current_finding'},
        {'assurance_scope': 'geometry_only'},
        {'context_policy': 'declared_frame'},
        {'padding': None},
        {'padding': 0},
        {'padding': '-1 nm'},
        {'padding': '1 fs'},
        {'view': 'other'},
        {'max_source_atoms': 5},
    ],
)
def test_unsupported_scope_or_invalid_request_refuses(source, options):
    params = {'method': 'enclosing_box', 'padding': '0.2 nm', **options}
    with pytest.raises((ArgumentError, CapabilityMismatchError)) as caught:
        declare().to_search_domain(source, **params)
    if 'method' in options:
        assert caught.value.extra['engine'] == 'bindingmt.site_domain_application'


@pytest.mark.parametrize('indices,padding', [([], '0.2 nm'), ([0], '0 nm')])
def test_empty_or_degenerate_support_is_not_repaired(source, indices, padding):
    with pytest.raises(ArgumentError):
        declare(indices).to_search_domain(
            source, method='enclosing_box', padding=padding
        )


@pytest.mark.parametrize(
    'failure',
    [
        MemoryError('authored'),
        OSError(errno.ENOMEM, 'authored'),
        KeyboardInterrupt(),
        SystemExit(7),
    ],
)
def test_coordinate_query_resource_errors_keep_original_identity(
    source, monkeypatch, failure
):
    get = msm.get

    def guarded(system, **kwargs):
        if kwargs.get('coordinates'):
            raise failure
        return get(system, **kwargs)

    monkeypatch.setattr(msm, 'get', guarded)
    with pytest.raises(type(failure)) as caught:
        domain(source)
    assert caught.value is failure


@pytest.mark.parametrize('field', ['center', 'lengths'])
def test_current_box_shape_cannot_pass_by_broadcasting(source, field):
    box = domain(source)
    getattr(box, field).shape = (1, 3)
    with pytest.raises(ArgumentError):
        box.to_dict()


@pytest.mark.parametrize(
    'path,value',
    [
        (('assurance', 'current_finding'), 'checked'),
        (('positions', 'sha256'), 'bad'),
        (('request', 'method'), 'sphere_enclosure'),
        (('consumer_binding', 'normalized_structure_index'), 2),
        (('losses', 'pbc'), True),
        (('record_id',), 'bad'),
        (('format_version',), '1.0'),
    ],
)
def test_resealed_invalid_application_records_refuse(source, path, value):
    record = domain(source).to_dict()
    node = record['derivation']
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    record['derivation']['content_sha256'] = _digest(record['derivation'])
    with pytest.raises((ArgumentError, CapabilityMismatchError)):
        bmt.BoxRegion.from_dict(record)


def test_consistent_reauthored_history_is_not_proof_at_use(source):
    record = domain(source).to_dict()
    recipe = record['derivation']
    for key in ['minimum', 'maximum']:
        old = recipe['bounds'][key]
        values = np.array(old['values'], dtype=float)
        values[0] += 1
        recipe['bounds'][key] = _sealed(values, key, fields=FIELDS)
    values = np.array(recipe['box']['center']['values'], dtype=float)
    values[0] += 1
    recipe['box']['center'] = _sealed(values, 'center', fields=FIELDS)
    record['center']['value'][0] += 1
    recipe['content_sha256'] = _digest(recipe)
    historical = bmt.BoxRegion.from_dict(record)
    with pytest.raises(ArgumentError):
        problem(source, historical)


@pytest.mark.parametrize('damage', ['missing', 'changed', 'boolean'])
def test_problem_reconstruction_requires_exact_construction_receipt(source, damage):
    record = problem(source, domain(source)).to_dict()
    if damage == 'missing':
        record.pop('search_domain_application')
    elif damage == 'boolean':
        record['search_domain_application']['chemical_state_index'] = True
    else:
        record['search_domain_application']['support_local_atom_indices'].reverse()
    with pytest.raises(ArgumentError):
        bmt.DockingProblem.from_dict(record, receptor=source, partner=MINIMAL_LIG_PDBQT)


def test_detached_receipt_and_empty_result_keep_historical_provenance(source, tmp_path):
    box = domain(source)
    p = problem(source, box)
    receipt = p.search_domain_application
    receipt['support_local_atom_indices'].clear()
    assert p.search_domain_application['support_local_atom_indices']
    result = bmt.DockingResult(
        [],
        problem_info=p.to_dict(),
        provenance={'scope': 'software_empty_result_control'},
    )
    record = result.to_dict()
    assert bmt.DockingResult.from_dict(record).to_dict() == record
    (tmp_path / 'application-result.json').write_text(
        json.dumps(record, sort_keys=True)
    )


def test_actual_public_vina_path_refuses_stale_input_before_preparation(
    source, monkeypatch
):
    import importlib.abc
    import sys

    from depdigest.core import checker

    from bindingmt.engines import vina as backend

    p = problem(source, domain(source))
    msm.set(
        p.receptor_molsys,
        selection=[1],
        structure_indices=[0],
        coordinates=puw.quantity([[[99.0, 0.0, 0.0]]], 'nm'),
    )

    def forbidden(*args, **kwargs):
        pytest.fail('Stale site application reached preparation or an engine')

    class ForbidEngine(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname.split('.')[0] in ('vina', 'meeko'):
                forbidden()

    monkeypatch.delitem(sys.modules, 'vina', raising=False)
    monkeypatch.setattr(sys, 'meta_path', [ForbidEngine(), *sys.meta_path])
    available = checker.is_installed
    monkeypatch.setattr(
        checker,
        'is_installed',
        lambda name: True if name == 'vina' else available(name),
    )
    monkeypatch.setattr(backend, 'prepare_receptor', forbidden)
    monkeypatch.setattr(backend, 'prepare_ligand', forbidden)
    with pytest.raises(ArgumentError):
        bmt.dock(p, bmt.VinaProtocol(cpu=1, seed=1))


def test_plain_box_keeps_old_export_and_problem_record(source):
    box = bmt.BoxRegion(center='[0,0,0] nm', lengths='[1,2,3] nm')
    p = problem(source, box)
    assert p.search_domain_application is None
    assert 'search_domain_application' not in p.to_dict()
    assert box.to_backend_box()['lengths'] == (10.0, 20.0, 30.0)


@pytest.mark.parametrize('unit', ['nm', 'pm', 'angstrom'])
def test_native_coordinate_units_preserve_actual_application(source, unit):
    coordinates = source.structures.coordinates + puw.quantity(
        [0.123456789, 0.234567891, 0.345678912], 'nm'
    )
    source.structures.coordinates = puw.convert(coordinates, to_unit=unit)
    box = domain(source)
    p = problem(source, box, receptor_selection=[4, 1, 3])
    np.testing.assert_allclose(
        box.to_backend_box(unit='nm', problem=p)['center'],
        [4.623456789, 1.484567891, -0.154321088],
    )


@pytest.mark.parametrize('producer,reader', [('pint', 'unyt'), ('unyt', 'pint')])
def test_decimal_application_recipe_has_faithful_cross_form_replay(
    source, producer, reader
):
    source.structures.coordinates = source.structures.coordinates + puw.quantity(
        [0.123456789, 0.234567891, 0.345678912], 'nm'
    )
    with puw.context(default_form=producer, standard_units=['nm', 'ps']):
        record = domain(source).to_dict()
    with puw.context(default_form=reader, standard_units=['pm', 'fs']):
        before = puw.configure.report()
        restored = bmt.BoxRegion.from_dict(record)
        replayed = restored.to_dict()
        assert replayed['derivation'] == record['derivation']
        for field in ['center', 'lengths']:
            np.testing.assert_array_equal(
                replayed[field]['value'], record[field]['value']
            )
            np.testing.assert_array_equal(
                puw.get_value(
                    puw.quantity(replayed[field]['value'], replayed[field]['unit']),
                    to_unit='nm',
                ),
                puw.get_value(
                    puw.quantity(record[field]['value'], record[field]['unit']),
                    to_unit='nm',
                ),
            )
        assert puw.configure.report() == before


def test_single_member_requires_explicit_positive_padding(source):
    box = domain(source, indices=[4])
    p = problem(source, box)
    np.testing.assert_allclose(
        box.to_backend_box(unit='nm', problem=p)['lengths'], [0.4, 0.4, 0.4]
    )


def prepared_control(p, retained=None):
    """An explicitly authored prepared view; no preparation method accuracy claim."""
    n_source = msm.get(p.receptor_molsys, n_atoms=True)
    retained = list(range(n_source)) if retained is None else retained
    n = len(retained)
    forward = [-1] * n_source
    for output, original in enumerate(retained):
        forward[original] = output
    return bmt.PreparedReceptor(
        state_id='authored-control',
        atom_names=[f'C{i}' for i in range(n)],
        group_names=['LIG'] * n,
        group_ids=[1] * n,
        coordinates=msm.get(p.receptor_molsys, selection=retained, coordinates=True)[0],
        atom_types=['C'] * n,
        charges=[0.0] * n,
        metadata={
            'molecular_projection': {
                'schema': 'molsysmt.atom_map@1',
                'n_source_atoms': n_source,
                'n_output_atoms': n,
                'source_to_output_atom_indices': forward,
                'output_to_source_atom_indices': retained,
            }
        },
    )


def test_real_owner_writer_and_staged_payload_are_admitted_before_engine_work(
    source, tmp_path
):
    from bindingmt.engines import vina as backend

    p = problem(source, domain(source))
    receptor = prepared_control(p)
    center, size = backend._vina_box(p.search_domain, problem=p)
    owned = []
    try:
        staged = object.__new__(bmt.VinaBackend)._stage_inputs(
            receptor,
            MINIMAL_LIG_PDBQT,
            bmt.VinaProtocol(cpu=1),
            owned,
            application_problem=p,
            application_box={'center': center, 'size': size, 'unit': 'angstrom'},
        )
        receipt = staged['artifacts']['receptor']['site_domain_application']
        assert receipt['support_pdbqt_atom_indices'] == [4, 1, 3]
        assert (
            receipt['receptor_payload_sha256']
            == staged['artifacts']['receptor']['sha256']
        )
        assert receipt['current_finding'] == 'not_assessed'
        (tmp_path / 'application-payload-use.json').write_text(
            json.dumps(receipt, sort_keys=True)
        )
    finally:
        backend._cleanup_staged_files(owned)


def test_prepared_coordinate_edits_refuse_with_actual_writer_bytes(source):
    from bindingmt.engines import vina as backend

    p = problem(source, domain(source))
    receptor = prepared_control(p)
    receptor.coordinates[1] += puw.quantity([0.1, 0.0, 0.0], 'nm')
    center, size = backend._vina_box(p.search_domain, problem=p)
    owned = []
    try:
        with pytest.raises(ArgumentError):
            object.__new__(bmt.VinaBackend)._stage_inputs(
                receptor,
                MINIMAL_LIG_PDBQT,
                bmt.VinaProtocol(cpu=1),
                owned,
                application_problem=p,
                application_box={'center': center, 'size': size, 'unit': 'angstrom'},
            )
    finally:
        backend._cleanup_staged_files(owned)


@pytest.mark.parametrize(
    'retained,accepted', [([0, 1, 2, 3, 4], True), ([0, 2, 3, 4, 5], False)]
)
def test_final_projection_requires_direct_retention_of_every_support_atom(
    source, retained, accepted
):
    from bindingmt.engines import vina as backend

    p = problem(source, domain(source))
    center, size = backend._vina_box(p.search_domain, problem=p)
    owned = []

    def stage():
        return object.__new__(bmt.VinaBackend)._stage_inputs(
            prepared_control(p, retained),
            MINIMAL_LIG_PDBQT,
            bmt.VinaProtocol(cpu=1),
            owned,
            application_problem=p,
            application_box={'center': center, 'size': size, 'unit': 'angstrom'},
        )

    try:
        if accepted:
            assert stage()['artifacts']['receptor']['site_domain_application'][
                'support_pdbqt_atom_indices'
            ] == [4, 1, 3]
        else:
            with pytest.raises(ArgumentError):
                stage()
    finally:
        backend._cleanup_staged_files(owned)


def test_actual_pdbqt_rounding_cannot_silently_escape_unpadded_box(source):
    from bindingmt.engines import vina as backend

    source.structures.coordinates = source.structures.coordinates + puw.quantity(
        [0.123456789, 0.234567891, 0.345678912], 'nm'
    )
    box = declare().to_search_domain(source, method='enclosing_box', padding='0 nm')
    p = problem(source, box)
    center, size = backend._vina_box(box, problem=p)
    owned = []
    try:
        with pytest.raises(ArgumentError) as caught:
            object.__new__(bmt.VinaBackend)._stage_inputs(
                prepared_control(p),
                MINIMAL_LIG_PDBQT,
                bmt.VinaProtocol(cpu=1),
                owned,
                application_problem=p,
                application_box={'center': center, 'size': size, 'unit': 'angstrom'},
            )
        assert caught.value.extra['arg_name'] == 'backend_box'
    finally:
        backend._cleanup_staged_files(owned)


def test_staged_bytes_must_match_the_actual_owner_export(source, monkeypatch):
    from bindingmt.engines import vina as backend

    p = problem(source, domain(source))
    center, size = backend._vina_box(p.search_domain, problem=p)
    stage = backend._stage_pdbqt_bytes
    monkeypatch.setattr(
        backend, '_stage_pdbqt_bytes', lambda data, files: stage(data + b'\n', files)
    )
    owned = []
    try:
        with pytest.raises(ArgumentError):
            object.__new__(bmt.VinaBackend)._stage_inputs(
                prepared_control(p),
                MINIMAL_LIG_PDBQT,
                bmt.VinaProtocol(cpu=1),
                owned,
                application_problem=p,
                application_box={'center': center, 'size': size, 'unit': 'angstrom'},
            )
    finally:
        backend._cleanup_staged_files(owned)


def test_replacing_a_valid_application_needs_a_new_problem(source):
    box = domain(source)
    p = problem(source, box)
    box._set_derivation(domain(source).derivation)
    with pytest.raises(ArgumentError):
        box.to_backend_box(problem=p)


@pytest.mark.parametrize('edit_prepared', [False, True])
def test_public_backend_reaches_engine_only_after_actual_staged_admission(
    source, monkeypatch, edit_prepared
):
    import sys
    from types import SimpleNamespace

    from depdigest.core import checker

    from bindingmt._private import site_application
    from bindingmt.engines import vina as backend

    p = problem(source, domain(source))
    receptor = prepared_control(p)
    if edit_prepared:
        receptor.coordinates[1] += puw.quantity([0.1, 0.0, 0.0], 'nm')
    admitted = []
    original = site_application.validate_staged_receptor

    def checked(*args, **kwargs):
        result = original(*args, **kwargs)
        admitted.append(result)
        return result

    class EngineReached(Exception):
        pass

    class SoftwareEngine:
        def __init__(self, **kwargs):
            assert admitted
            raise EngineReached('software sentinel; native engine is never constructed')

    monkeypatch.setattr(site_application, 'validate_staged_receptor', checked)
    monkeypatch.setattr(backend, 'prepare_receptor', lambda *args, **kwargs: receptor)
    available = checker.is_installed
    monkeypatch.setattr(
        checker,
        'is_installed',
        lambda name: True if name == 'vina' else available(name),
    )
    monkeypatch.setitem(
        sys.modules,
        'vina',
        SimpleNamespace(Vina=SoftwareEngine, __version__='software-control'),
    )
    with pytest.raises(ArgumentError if edit_prepared else EngineReached):
        bmt.dock(p, bmt.VinaProtocol(cpu=1, seed=1, allow_provisional_preparation=True))
    assert bool(admitted) is not edit_prepared
