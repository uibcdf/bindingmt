"""Frozen selection snapshots keep owner semantics and explicit chemical scope."""

import errno
import importlib
import json
from copy import deepcopy

import molsysmt as msm
import numpy as np
import pandas as pd
import pytest
from molsysmt._private.smonitor import ArgumentError as OwnerArgumentError
from molsysmt.native import MolSys

import bindingmt as bmt
from bindingmt._private.interface_records import _digest
from bindingmt._private.smonitor import ArgumentError, CapabilityMismatchError

puw = msm.pyunitwizard
SOURCE = {'owner': 'authored.control', 'id': 'selection-source', 'revision': '1'}


@pytest.fixture
def source():
    system = MolSys(n_atoms=4, n_groups=3, n_chains=1, n_molecules=3, n_entities=3)
    top = system.topology
    top.atoms['atom_id'] = ['02', 'same', '02', 'last']
    top.atoms['atom_name'] = ['a', 'b', 'c', 'd']
    top.atoms['atom_type'] = ['C'] * 4
    top.atoms['group_index'] = pd.array([0, 1, 0, 1], dtype='Int64')
    top.atoms['chain_index'] = pd.array([0] * 4, dtype='Int64')
    top.groups['group_id'] = ['02', 'same', '02']
    top.groups['group_name'] = ['LIG'] * 3
    top.groups['group_type'] = ['small molecule'] * 3
    top.groups['molecule_index'] = pd.array([0, 1, 2], dtype='Int64')
    top.molecules['entity_index'] = pd.array([0, 1, 2], dtype='Int64')
    top.chains['chain_id'] = ['A']
    coordinates = np.array(
        [
            [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [0.2, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        ]
    )
    system.structures.coordinates = puw.quantity(coordinates, 'nm')
    system.structures.time = puw.quantity([0.0, 1.0, 2.0], 'ps')
    system.chemical_states.append_state()
    msm.set(system, chemical_state=0, formal_charge=[1, 0, 0, 0])
    msm.set(system, chemical_state=1, formal_charge=[0, 1, 0, 0])
    msm.set(system, structure_chemical_state_index=[0, 1, 0])
    return system


def test_owner_resident_inventory_composition_needs_no_structures(source):
    inventory = [source.topology, source.chemical_states]
    assert msm.select(
        inventory,
        'group_index in [1, 0]',
        element='group',
        structure_indices=[],
        chemical_state=0,
    ) == [0, 1]
    assert msm.select(
        inventory, 'formal_charge>0', structure_indices=[], chemical_state=0
    ) == [0]
    assert msm.select(
        inventory, 'formal_charge>0', structure_indices=[], chemical_state=1
    ) == [1]
    assert msm.get(source, reference_chemical_state_index=True) == 0


@pytest.mark.parametrize('structure,expected', [(0, [0, 1]), (1, [0, 2]), (2, [0])])
def test_owner_spatial_query_respects_one_explicit_structure(
    source, structure, expected
):
    assert (
        msm.select(
            source,
            'all within 0.3 nm of atom_index==0',
            structure_indices=[structure],
            chemical_state='structure',
        )
        == expected
    )


def make(source, selection='all', **options):
    return bmt.BindingSite.from_selection(source, selection, source=SOURCE, **options)


@pytest.mark.parametrize('kind', ['object', 'path', 'mapping', 'domains'])
def test_native_admission_selection_refuses_foreign_inputs_without_discovery(
    source, monkeypatch, kind
):
    foreign = {
        'object': object(),
        'path': 'not-loaded.h5msm',
        'mapping': {'coordinates': []},
        'domains': [source.topology, source.chemical_states],
    }[kind]

    def forbidden(*args, **kwargs):
        pytest.fail('Native-only selection attempted foreign form discovery or queries')

    for name in ['get_form', 'get', 'select', 'convert', 'extract']:
        monkeypatch.setattr(msm, name, forbidden)
    with pytest.raises(CapabilityMismatchError) as caught:
        make(foreign)
    assert caught.value.extra['capability'] == 'resident native selection source'


@pytest.mark.parametrize(
    'query,expected',
    [
        ('all', [0, 1, 2, 3]),
        ([3, 1, 0], [3, 1, 0]),
        ((2, 0), [2, 0]),
        (np.array([3, 0], dtype=np.int32), [3, 0]),
        (np.int64(2), [2]),
        ([], []),
        ('atom_id=="02"', [0, 2]),
        ('atom_name=="absent"', []),
    ],
)
def test_fixed_inventory_snapshot_keeps_request_order_and_original_source(
    source, query, expected
):
    before = source.topology.atoms.copy(deep=True)
    coordinates = puw.get_value(source.structures.coordinates, to_unit='nm').copy()
    declaration = make(source, query)
    record = declaration.to_dict()
    assert record['format_version'] == '0.4'
    view = record['views'][0]
    assert view['payload']['source_indices'] == expected
    assert view['role'] == 'selected_members'
    assert view['payload']['resolution'] == {
        'state': 'resolved',
        'cause': 'captured_selection',
    }
    capture = view['payload']['selection']
    assert capture['effective']['mode'] == 'inventory'
    assert capture['effective']['structure_indices'] == []
    assert capture['effective']['chemical_state_index'] == 0
    assert capture['effective']['chemical_state_policy'] == 'reference'
    assert capture['membership_policy'] == 'fixed_snapshot'
    assert capture['current_applicability'] == 'not_assessed'
    assert declaration.definition['scope']['kind'] == 'independent'
    assert declaration.definition['binding_support'] == 'unknown'
    assert declaration.inspect()['selection_capture'] == capture
    assert declaration.inspect()['molecular_resolution'] == 'not_performed'
    if not expected:
        assert (
            declaration.resolve_region(source, structure_index=0)['projection'][
                'coverage'
            ]
            == 'declared_empty'
        )
    assert bmt.BindingSite.from_dict(record, experimental=True).to_dict() == record
    assert declaration.copy().to_dict() == record
    pd.testing.assert_frame_equal(source.topology.atoms, before)
    np.testing.assert_array_equal(
        puw.get_value(source.structures.coordinates, to_unit='nm'), coordinates
    )


def test_inventory_queries_use_no_structures_or_conversion(source, monkeypatch):
    select = msm.select
    seen = []

    def guarded(system, **kwargs):
        assert system == [source.topology, source.chemical_states]
        assert kwargs['structure_indices'] == []
        assert kwargs['chemical_state'] == 0
        seen.append(kwargs['selection'])
        return select(system, **kwargs)

    def forbidden(*args, **kwargs):
        pytest.fail('Inventory capture requested conversion or coordinates')

    monkeypatch.setattr(msm, 'select', guarded)
    monkeypatch.setattr(msm, 'convert', forbidden)
    monkeypatch.setattr(msm.structure, 'get_contacts', forbidden)
    declaration = make(source, 'formal_charge>0')
    assert declaration.views[0]['payload']['source_indices'] == [0]
    assert seen == ['formal_charge>0', 'formal_charge>0']


def test_group_inventory_preserves_known_empty_group_without_expanding_it(source):
    declaration = make(source, 'all', element='group')
    assert declaration.views[0]['payload']['source_indices'] == [0, 1, 2]
    region = declaration.resolve_region(
        source, structure_index=0, expansion='source_members'
    )
    assert [
        entry['target_atom_indices'] for entry in region['projection']['entries']
    ] == [[0, 2], [1, 3], []]
    assert region['projection']['entries'][2]['state'] == 'empty'


@pytest.mark.parametrize('state,expected', [(0, [0]), (1, [1])])
def test_inventory_explicit_chemical_dependency_is_preserved_in_projection(
    source, state, expected
):
    declaration = make(source, 'formal_charge>0', chemical_state=state)
    assert declaration.views[0]['payload']['source_indices'] == expected
    region = declaration.resolve_region(source, structure_index=state)
    assert region['projection']['exported_target_atom_indices'] == expected
    with pytest.raises(ArgumentError):
        declaration.resolve_region(source, structure_index=1 - state)


@pytest.mark.parametrize('structure,expected', [(0, [0]), (1, [1]), (2, [0])])
def test_structure_resolved_selection_records_actual_state_and_scope(
    source, structure, expected
):
    declaration = make(source, 'formal_charge>0', structure_indices=[structure])
    payload = declaration.views[0]['payload']
    assert payload['source_indices'] == expected
    assert payload['selection']['effective']['chemical_state_policy'] == 'structure'
    assert payload['selection']['effective']['chemical_state_index'] == (
        1 if structure == 1 else 0
    )
    assert declaration.definition['scope']['source_structure_indices'] == [structure]
    region = declaration.resolve_region(source)
    assert region['projection']['exported_target_atom_indices'] == expected


@pytest.mark.parametrize('structure,expected', [(0, [0, 1]), (1, [0, 2]), (2, [0])])
@pytest.mark.parametrize(
    'query',
    ['all within 0.3 nm of atom_index==0', 'all within 3 angstrom of atom_index==0'],
)
def test_spatial_snapshot_uses_owner_queries_and_retains_physical_expression(
    source, structure, expected, query
):
    declaration = make(source, query, structure_indices=structure)
    payload = declaration.views[0]['payload']
    assert payload['source_indices'] == expected
    assert payload['selection']['request']['selection'] == query
    assert payload['selection']['effective']['budgeted_dense_pairs'] == 16
    assert (
        declaration.resolve_region(source)['projection']['exported_target_atom_indices']
        == expected
    )


def test_spatial_query_without_structure_raises_owner_diagnostic_instead_of_empty(
    source,
):
    with pytest.raises(OwnerArgumentError) as caught:
        make(source, 'all within 0.3 nm of atom_index==0')
    assert 'coordinates' in str(caught.value).lower()


def test_selection_is_a_fixed_snapshot_after_coordinate_edits(source, monkeypatch):
    declaration = make(
        source, 'all within 0.3 nm of atom_index==0', structure_indices=[0]
    )
    before = declaration.to_dict()
    msm.set(
        source,
        selection=[1],
        structure_indices=[0],
        coordinates=puw.quantity([[[5.0, 0.0, 0.0]]], 'nm'),
    )
    assert msm.select(
        source,
        'all within 0.3 nm of atom_index==0',
        structure_indices=[0],
        chemical_state='structure',
    ) == [0]
    original = msm.select

    def guarded(system, selection='all', **kwargs):
        assert not isinstance(selection, str), (
            'Projection reevaluated a historical predicate'
        )
        return original(system, selection=selection, **kwargs)

    monkeypatch.setattr(msm, 'select', guarded)
    assert declaration.resolve_region(source)['projection'][
        'exported_target_atom_indices'
    ] == [0, 1]
    assert declaration.to_dict() == before


def test_capture_survives_immutable_claim_revision_and_saved_result_replay(
    source, tmp_path
):
    declaration = make(source, [3, 1], name='candidate')
    original = declaration.to_dict()
    revised = declaration.with_definition(
        claim='New interpretation', criterion={'kind': 'authored'}, reason='Review'
    )
    assert (
        revised.views == declaration.views
        and revised.record_id != declaration.record_id
    )
    assert declaration.to_dict() == original
    path = tmp_path / 'selected-site.json'
    declaration.save(path)
    restored = bmt.BindingSite.read(path, experimental=True)
    assert restored.to_dict() == original
    region = declaration.resolve_region(source, structure_index=0)
    (tmp_path / 'selected-region.json').write_text(json.dumps(region, sort_keys=True))
    assert bmt.BindingSite.region_from_dict(region, experimental=True) == region


def test_historical_site_and_region_readers_never_rerun_selection(source, monkeypatch):
    declaration = make(source, 'formal_charge>0', structure_indices=[0])
    region = declaration.resolve_region(source)

    def forbidden(*args, **kwargs):
        pytest.fail('Historical reader accessed a molecular source')

    for name in ['select', 'get', 'get_form', 'convert', 'extract']:
        monkeypatch.setattr(msm, name, forbidden)
    assert (
        bmt.BindingSite.from_dict(declaration.to_dict(), experimental=True).to_dict()
        == declaration.to_dict()
    )
    assert bmt.BindingSite.region_from_dict(region, experimental=True) == region


@pytest.mark.parametrize(
    'options',
    [
        {'element': 'residue'},
        {'syntax': 'MDTraj'},
        {'role': 'declared_members'},
        {'structure_indices': []},
        {'structure_indices': [0, 1]},
        {'structure_indices': [2, 2]},
        {'structure_indices': [-1]},
        {'structure_indices': [3]},
        {'structure_indices': True},
        {'structure_indices': 1, 'chemical_state': 'reference'},
        {'structure_indices': 1, 'chemical_state': 0},
        {'chemical_state': 'structure'},
        {'chemical_state': -1},
        {'chemical_state': True},
        {'chemical_state': 2},
        {'chemical_state': 'other'},
        {'max_source_atoms': 3},
        {'max_structure_occurrences': 2},
        {'structure_indices': 0, 'max_dense_pairs': 15},
    ],
)
def test_wrong_context_or_budget_refuses(source, options):
    with pytest.raises((ArgumentError, CapabilityMismatchError)):
        make(source, 'formal_charge>0', **options)


@pytest.mark.parametrize(
    'query',
    [
        None,
        True,
        1.0,
        -1,
        [True],
        [0, 0],
        [[0]],
        ['all'],
        np.array([[0, 1]]),
        'atom_index in @outside',
    ],
)
def test_nonportable_or_ambiguous_selections_refuse(source, query):
    with pytest.raises((ArgumentError, CapabilityMismatchError)):
        make(source, query)


@pytest.mark.parametrize(
    'failure',
    [
        MemoryError('owned'),
        OSError(errno.ENOMEM, 'owned'),
        KeyboardInterrupt(),
        SystemExit(3),
    ],
)
def test_direct_owner_query_errors_keep_original_identity(source, monkeypatch, failure):
    def failed(*args, **kwargs):
        raise failure

    monkeypatch.setattr(msm, 'select', failed)
    with pytest.raises(type(failure)) as caught:
        make(source)
    assert caught.value is failure


@pytest.mark.parametrize('mode', ['inventory', 'single_structure'])
@pytest.mark.parametrize('element', ['atom', 'group'])
@pytest.mark.parametrize('failed_call', [1, 2])
@pytest.mark.parametrize('kind', ['memory', 'enomem', 'interrupt', 'exit'])
def test_actual_expression_boundary_preserves_original_failure(
    source, monkeypatch, mode, element, failed_call, kind
):
    # Instrument only the owner's evaluator, leaving both public APIs intact.
    # Inject exceptions instead of inducing allocation pressure or real signals.
    selector = importlib.import_module('molsysmt.basic.selector')
    original = selector._dict_select['MolSysMT']
    failure = {
        'memory': MemoryError('synthetic expression capacity'),
        'enomem': OSError(errno.ENOMEM, 'synthetic expression capacity'),
        'interrupt': KeyboardInterrupt('synthetic expression interrupt'),
        'exit': SystemExit(19),
    }[kind]
    cause = RuntimeError('original expression cause')
    failure.__cause__ = cause
    failure.add_note('Original expression note')
    calls = 0

    def injected(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == failed_call:
            raise failure
        return original(*args, **kwargs)

    monkeypatch.setitem(selector._dict_select, 'MolSysMT', injected)
    options = (
        {'structure_indices': 1, 'chemical_state': 'structure'}
        if mode == 'single_structure'
        else {'chemical_state': 1}
    )
    with pytest.raises(type(failure)) as caught:
        make(source, 'atom_index>=0', element=element, **options)
    assert calls == failed_call
    assert caught.value is failure
    assert caught.value.__cause__ is cause
    assert caught.value.__notes__ == ['Original expression note']


@pytest.mark.parametrize('mode', ['inventory', 'single_structure'])
@pytest.mark.parametrize('kind', ['memory', 'enomem'])
def test_actual_pandas_evaluation_resource_failure_reaches_public_site(
    source, monkeypatch, mode, kind
):
    failure = (
        MemoryError('synthetic dataframe capacity')
        if kind == 'memory'
        else OSError(errno.ENOMEM, 'synthetic dataframe capacity')
    )
    calls = []

    def injected(frame, expression, *args, **kwargs):
        calls.append(expression)
        raise failure

    monkeypatch.setattr(pd.DataFrame, 'query', injected)
    options = {'structure_indices': 1} if mode == 'single_structure' else {}
    with pytest.raises(type(failure)) as caught:
        make(source, 'atom_name=="a"', **options)
    assert len(calls) == 1
    assert caught.value is failure
    assert caught.value.__cause__ is None


@pytest.mark.parametrize('mode', ['inventory', 'single_structure'])
@pytest.mark.parametrize('element', ['atom', 'group'])
@pytest.mark.parametrize('kind', ['value', 'syntax', 'other_oserror'])
def test_actual_expression_parse_failures_retain_owner_diagnostic(
    source, monkeypatch, mode, element, kind
):
    selector = importlib.import_module('molsysmt.basic.selector')
    failure = {
        'value': ValueError('synthetic invalid expression'),
        'syntax': SyntaxError('synthetic invalid expression'),
        'other_oserror': OSError(errno.ENOENT, 'synthetic parser input missing'),
    }[kind]
    calls = []

    def injected(*args, **kwargs):
        calls.append(args)
        raise failure

    monkeypatch.setitem(selector._dict_select, 'MolSysMT', injected)
    options = {'structure_indices': 1} if mode == 'single_structure' else {}
    with pytest.raises(OwnerArgumentError) as caught:
        make(source, 'atom_index>=0', element=element, **options)
    assert len(calls) == 1
    assert caught.value.__cause__ is failure


def test_changed_result_during_capture_refuses(source, monkeypatch):
    original = msm.select
    calls = 0

    def changing(*args, **kwargs):
        nonlocal calls
        calls += 1
        result = original(*args, **kwargs)
        return result if calls == 1 else list(reversed(result))

    monkeypatch.setattr(msm, 'select', changing)
    with pytest.raises(ArgumentError):
        make(source)


@pytest.mark.parametrize(
    'path,value',
    [
        (('selection', 'effective', 'chemical_state_index'), 1),
        (('selection', 'effective', 'mode'), 'single_structure'),
        (('selection', 'observations', 'n_atoms'), 3),
        (('selection', 'request', 'selection'), [1, 2]),
        (('selection', 'request', 'chemical_state'), 'structure'),
        (('selection', 'provider', 'library'), 'bindingmt'),
        (('selection', 'membership_policy'), 'live_rule'),
        (('selection', 'current_applicability'), 'verified'),
        (('resolution', 'state'), 'declared'),
        (('source_indices',), [3.0, 1]),
        (('index_space',), 'target_atom_index'),
    ],
)
def test_resealed_capture_inconsistencies_refuse(source, path, value):
    record = make(source, [3, 1]).to_dict()
    node = record['views'][0]['payload']
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = deepcopy(value)
    record['content_sha256'] = _digest(record)
    with pytest.raises((ArgumentError, CapabilityMismatchError)):
        bmt.BindingSite.from_dict(record, experimental=True)


@pytest.mark.parametrize('form', ['pint', 'unyt'])
def test_unit_form_and_application_policy_are_preserved(source, form):
    with puw.context(default_form=form, standard_units=['pm', 'fs']):
        before = puw.configure.report()
        declaration = make(
            source, 'all within 3 angstrom of atom_index==0', structure_indices=0
        )
        assert declaration.views[0]['payload']['source_indices'] == [0, 1]
        assert puw.configure.report() == before


def test_known_captured_axis_changes_refuse_before_projection(source):
    declaration = make(source)
    smaller = msm.extract(source, selection=[0, 1, 2])
    with pytest.raises(ArgumentError):
        declaration.resolve_region(smaller, structure_index=0)


def test_inventory_selection_works_without_any_coordinate_structure(source):
    inventory_source = msm.extract(source, structure_indices=[])
    assert msm.get(inventory_source, n_structures=True) == 0
    declaration = make(inventory_source, 'formal_charge>0')
    assert declaration.views[0]['payload']['source_indices'] == [0]
    assert (
        declaration.views[0]['payload']['selection']['effective']['mode'] == 'inventory'
    )


def test_unknown_structure_association_never_falls_back_to_reference(source):
    msm.set(source, structure_indices=[0], structure_chemical_state_index=pd.NA)
    with pytest.raises(ArgumentError):
        make(source, 'formal_charge>0', structure_indices=[0])


def test_pair_budget_is_checked_before_any_owner_selection(source, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Over-budget request reached owner selection')

    monkeypatch.setattr(msm, 'select', forbidden)
    with pytest.raises(ArgumentError):
        make(
            source,
            'all within 0.3 nm of atom_index==0',
            structure_indices=[0],
            max_dense_pairs=15,
        )


def test_resealed_structure_scope_cannot_diverge_from_capture(source):
    record = make(source, 'formal_charge>0', structure_indices=[0]).to_dict()
    scope = record['views'][0]['scope']
    scope['source_structure_indices'] = [2]
    record['definition']['scope'] = deepcopy(scope)
    record['content_sha256'] = _digest(record)
    with pytest.raises(ArgumentError):
        bmt.BindingSite.from_dict(record, experimental=True)
