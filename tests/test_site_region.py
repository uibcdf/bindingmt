"""Native site projections preserve declared axes without inventing identity."""

import errno
import json
from copy import deepcopy
from pathlib import Path

import molsysmt as msm
import numpy as np
import pandas as pd
import pytest
from molsysmt._private.smonitor import StructuralInconsistencyError
from molsysmt.native import MolSys

import bindingmt as bmt
from bindingmt._private.interface_records import _digest
from bindingmt._private.smonitor import ArgumentError, CapabilityMismatchError

puw = msm.pyunitwizard
SOURCE = {'owner': 'authored.control', 'id': 'source-A', 'revision': '1'}


@pytest.fixture
def source():
    """A transport control, not a prepared or scientifically qualified case."""
    molsys = MolSys(n_atoms=6, n_groups=3, n_chains=1, n_molecules=3, n_entities=3)
    top = molsys.topology
    top.atoms['atom_id'] = ['02', 'same', '02', '4', '05', 'last']
    top.atoms['atom_name'] = list('abcdef')
    top.atoms['atom_type'] = ['C'] * 6
    top.atoms['group_index'] = pd.array([0, 1, 0, 2, 1, 2], dtype='Int64')
    top.atoms['chain_index'] = pd.array([0] * 6, dtype='Int64')
    top.groups['group_id'] = ['02', 'same', '02']
    top.groups['group_name'] = ['LIG'] * 3
    top.groups['group_type'] = ['small molecule'] * 3
    top.groups['molecule_index'] = pd.array([0, 1, 2], dtype='Int64')
    top.molecules['entity_index'] = pd.array([0, 1, 2], dtype='Int64')
    top.chains['chain_id'] = ['A']
    molsys.structures.coordinates = puw.quantity(np.arange(54.0).reshape(3, 6, 3), 'nm')
    molsys.structures.time = puw.quantity([10.0, 20.0, 30.0], 'ps')
    molsys.chemical_states.append_state()
    msm.set(molsys, chemical_state=0, formal_charge=[1, 0, 0, 0, 0, 0])
    msm.set(molsys, chemical_state=1, formal_charge=[0, 0, 1, 0, 0, 0])
    msm.set(molsys, structure_chemical_state_index=[1, 0, 0])
    return molsys


def site(element='atom', indices=(5, 1, 2), **options):
    return getattr(
        bmt.BindingSite, f'from_{"atoms" if element == "atom" else "groups"}'
    )(source=SOURCE, **{f'{element}_indices': list(indices)}, **options)


def derived(source):
    subset = msm.extract(source, selection=[5, 1, 2], structure_indices=[2, 0, 2])
    first = dict(
        schema='molsysmt.atom_map@1',
        n_source_atoms=6,
        n_output_atoms=3,
        source_to_output_atom_indices=[-1, 0, 1, -1, -1, 2],
        output_to_source_atom_indices=[1, 2, 5],
    )
    result = msm.permute_atoms(subset, [2, 0, 1], return_report=True)
    return dict(
        target=result['molecular_system'],
        atom_map=msm.compose_atom_maps([first, result['report']]),
        structure_correspondence=[2, 0, 2],
        structure_index=2,
        target_structure_index=0,
    )


def test_direct_atom_order_detached_record_and_unchanged_declaration(source):
    declaration = site(structure_indices=[2, 0, 2], chemical_state_indices=[0, 1, 0])
    before = declaration.to_dict()
    coordinates = puw.get_value(source.structures.coordinates, to_unit='nm').copy()
    result = declaration.resolve_region(source, structure_index=2)
    projection = result['projection']
    assert projection['exported_target_atom_indices'] == [5, 1, 2]
    assert [entry['source_index'] for entry in projection['entries']] == [5, 1, 2]
    assert projection['axes']['declaration_occurrences'] == [0, 2]
    assert projection['axes']['chemical_state_index'] == 0
    assert msm.select(
        source, 'formal_charge>0', structure_indices=[2], chemical_state='structure'
    ) == [0]
    assert 'binding support' in projection['assurance']['not_assessed']
    assert declaration.to_dict() == before
    assert declaration.inspect()['molecular_resolution'] == 'not_performed'
    assert declaration.inspect()['capabilities']['molecular_projection'] is True
    np.testing.assert_array_equal(
        puw.get_value(source.structures.coordinates, to_unit='nm'), coordinates
    )
    roundtrip = bmt.BindingSite.region_from_dict(
        json.loads(json.dumps(result, sort_keys=True)), experimental=True
    )
    assert roundtrip == result
    roundtrip['projection']['entries'].clear()
    assert result['projection']['entries']


@pytest.mark.parametrize('expansion', ['source_members', 'target_group_members'])
def test_direct_group_expansion_is_entry_aligned(source, expansion):
    result = site('group', [2, 0]).resolve_region(
        source, structure_index=2, expansion=expansion
    )
    assert result['observations']['source']['group_members'] == [[3, 5], [0, 2]]
    assert result['projection']['exported_target_atom_indices'] == [3, 5, 0, 2]
    assert bmt.BindingSite.region_from_dict(result, experimental=True) == result


def test_derived_permutation_repeated_structures_and_state_positions(source):
    options = derived(source)
    declaration = site(structure_indices=[2, 0, 2], chemical_state_indices=[0, 1, 0])
    result = declaration.resolve_region(source, **options)
    p = result['projection']
    assert p['exported_target_atom_indices'] == [0, 1, 2]
    assert [entry['output_source_atom_indices'] for entry in p['entries']] == [
        [5],
        [1],
        [2],
    ]
    assert p['axes']['target_to_source'] == [2, 0, 2]
    assert p['axes']['source_to_target'] == [[1], [], [0, 2]]
    assert (
        p['axes']['source_structure_index'] == 2
        and p['axes']['target_structure_index'] == 0
    )
    assert p['axes']['chemical_state_index'] == 0
    assert p['group_correspondence'] == [1, 0, 2]
    assert bmt.BindingSite.region_from_dict(result, experimental=True) == result
    options.pop('target_structure_index')
    with pytest.raises(ArgumentError):
        declaration.resolve_region(source, **options)
    other = declaration.resolve_region(source, **options, target_structure_index=2)
    assert other['projection']['axes']['target_structure_index'] == 2


@pytest.mark.parametrize('expansion', ['source_members', 'target_group_members'])
def test_known_omissions_are_retained_in_partial_group_results(source, expansion):
    options = derived(source)
    declaration = site('group', [0, 2])
    with pytest.raises(ArgumentError):
        declaration.resolve_region(source, **options, expansion=expansion)
    result = declaration.resolve_region(
        source, **options, expansion=expansion, coverage_policy='allow_partial'
    )
    assert [
        entry['missing_source_atom_indices']
        for entry in result['projection']['entries']
    ] == [[0], [3]]
    assert [
        entry['target_atom_indices'] for entry in result['projection']['entries']
    ] == [[2], [0]]
    assert result['projection']['coverage'] == 'partial'
    assert bmt.BindingSite.region_from_dict(result, experimental=True) == result


def test_partial_atom_entries_never_pass_minus_one_to_a_source_lookup(source):
    options = derived(source)
    declaration = site(indices=[0, 5])
    with pytest.raises(ArgumentError):
        declaration.resolve_region(source, **options)
    result = declaration.resolve_region(
        source, **options, coverage_policy='allow_partial'
    )
    first, last = result['projection']['entries']
    assert (
        first['missing_source_atom_indices'] == [0]
        and first['target_atom_indices'] == []
    )
    assert last['target_atom_indices'] == [0]
    assert result['projection']['coverage'] == 'partial'


def test_generated_group_atom_never_inherits_original_identity():
    builder = msm.MolSysBuilder()
    builder.add_atom(atom_type='C', atom_name='parent')
    builder.set_coordinates(puw.quantity([[[0.0, 0.0, 0.0]]], 'nm'))
    source = builder.build()
    addition = msm.build.add_terminal_atoms(
        source,
        [{'parent_atom_index': 0, 'atom_type': 'H'}],
        puw.quantity([[[0.0, 0.0, 0.109]]], 'nm'),
    )
    options = dict(
        target=addition['molecular_system'],
        atom_map=addition['report'],
        structure_correspondence=[0],
    )
    declaration = site('group', [0])
    original = declaration.resolve_region(source, **options, expansion='source_members')
    current = declaration.resolve_region(
        source, **options, expansion='target_group_members'
    )
    assert original['projection']['exported_target_atom_indices'] == [0]
    assert current['projection']['exported_target_atom_indices'] == [0, 1]
    entry = current['projection']['entries'][0]
    assert entry['output_source_atom_indices'] == [0, None]
    assert entry['generated_target_atom_indices'] == [1]
    assert current['site'] == declaration.to_dict()
    assert bmt.BindingSite.region_from_dict(current, experimental=True) == current


@pytest.mark.parametrize('element', ['atom', 'group'])
def test_declared_empty_does_not_become_a_negative_finding(source, element):
    result = site(element, []).resolve_region(
        source,
        structure_index=2,
        expansion=None if element == 'atom' else 'source_members',
    )
    assert result['projection']['entries'] == []
    assert result['projection']['coverage'] == 'declared_empty'
    assert result['site']['definition']['binding_support'] == 'unknown'


def test_known_empty_group_and_unknown_parentage_are_distinct(source):
    source.topology.atoms['group_index'] = pd.array([0, 1, 0, 1, 1, 1], dtype='Int64')
    result = site('group', [2]).resolve_region(
        source, structure_index=2, expansion='source_members'
    )
    assert result['projection']['entries'][0]['state'] == 'empty'
    source.topology.atoms.loc[5, 'group_index'] = pd.NA
    with pytest.raises(CapabilityMismatchError):
        site('group', [2]).resolve_region(
            source,
            structure_index=2,
            expansion='source_members',
            coverage_policy='allow_partial',
        )


@pytest.mark.parametrize(
    'failure',
    [
        MemoryError('owned'),
        OSError(errno.ENOMEM, 'owned'),
        KeyboardInterrupt(),
        SystemExit(7),
    ],
)
def test_provider_resource_and_process_errors_keep_original_identity(
    source, monkeypatch, failure
):
    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(msm, 'get', fail)
    with pytest.raises(type(failure)) as caught:
        site().resolve_region(source, structure_index=2)
    assert caught.value is failure


def test_frozen_membership_record_is_consumed_without_rewriting_its_producer(source):
    path = Path(__file__).parent / 'data/site_geometry_compat/e1e6a33/atoms.json'
    before = path.read_bytes()
    declaration = bmt.BindingSite.from_dict(json.loads(before), experimental=True)
    result = declaration.resolve_region(source, structure_index=2)
    assert result['projection']['exported_target_atom_indices'] == [5, 1, 3]
    assert result['site']['producer'] == json.loads(before)['producer']
    assert path.read_bytes() == before


@pytest.mark.parametrize(
    'path,value',
    [
        (('projection', 'exported_target_atom_indices'), [1, 5, 2]),
        (('projection', 'exported_target_atom_indices'), [5.0, 1, 2]),
        (('projection', 'coverage'), 'partial'),
        (('projection', 'axes', 'chemical_state_index'), 1),
        (('definition_ref', 'record_id'), 'wrong'),
        (('view_ref', 'key'), 'wrong'),
        (('observations', 'target', 'atom_group_indices'), [0, 1, 0, 2, 1, 0]),
    ],
)
def test_resealed_projection_changes_refuse(source, path, value):
    record = site().resolve_region(source, structure_index=2)
    target = record
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    record['content_sha256'] = _digest(record)
    with pytest.raises(ArgumentError):
        bmt.BindingSite.region_from_dict(record, experimental=True)


def test_replay_is_detached_and_never_queries_a_molecular_source(
    source, monkeypatch, tmp_path
):
    result = site().resolve_region(source, structure_index=2)
    path = tmp_path / 'region.json'
    path.write_text(json.dumps(result, sort_keys=True))

    def forbidden(*args, **kwargs):
        pytest.fail('Historical replay queried a molecular source')

    for name in ['get', 'select', 'extract', 'convert', 'permute_atoms', 'get_form']:
        monkeypatch.setattr(msm, name, forbidden)
    assert (
        bmt.BindingSite.region_from_dict(
            json.loads(path.read_text()), experimental=True
        )
        == result
    )
    with pytest.raises(ArgumentError):
        bmt.BindingSite.region_from_dict(result)


@pytest.mark.parametrize(
    'options',
    [
        {'structure_index': -1},
        {'structure_index': True},
        {'structure_index': 3},
        {'structure_index': 2, 'chemical_state': 1},
        {'structure_index': 2, 'chemical_state': 'reference'},
        {'structure_index': 2, 'context_policy': 'declared_frame'},
        {'structure_index': 2, 'coverage_policy': 'ignore'},
        {'structure_index': 2, 'view': 'absent'},
        {'structure_index': 2, 'expansion': 'all_atoms'},
        {'structure_index': 2, 'max_source_atoms': 5},
        {'structure_index': 2, 'max_structure_occurrences': 2},
    ],
)
def test_unsupported_or_invalid_choices_refuse(source, options):
    with pytest.raises((ArgumentError, CapabilityMismatchError)):
        site().resolve_region(source, **options)


def test_required_structure_scope_and_group_expansion_are_explicit(source):
    with pytest.raises(ArgumentError):
        site().resolve_region(source)
    with pytest.raises(CapabilityMismatchError):
        site('group', [2, 0]).resolve_region(source, structure_index=2)
    declaration = site(structure_indices=[2, 2])
    assert (
        declaration.resolve_region(source)['projection']['axes'][
            'source_structure_index'
        ]
        == 2
    )


def test_explicit_target_never_uses_implicit_identity(source):
    options = derived(source)
    for missing in ('atom_map', 'structure_correspondence'):
        incomplete = {k: v for k, v in options.items() if k != missing}
        with pytest.raises(ArgumentError):
            site().resolve_region(source, **incomplete)


def test_selected_unknown_state_refuses_without_reference_fallback(source):
    msm.set(source, structure_indices=[2], structure_chemical_state_index=pd.NA)
    with pytest.raises(StructuralInconsistencyError):
        site().resolve_region(source, structure_index=2)


def test_observed_change_during_check_refuses(source, monkeypatch):
    original = msm.get
    calls = 0

    def changing(system, **kwargs):
        nonlocal calls
        if kwargs.get('group_index'):
            calls += 1
            if calls == 3:
                source.topology.atoms.loc[0, 'group_index'] = 1
        return original(system, **kwargs)

    monkeypatch.setattr(msm, 'get', changing)
    with pytest.raises(ArgumentError):
        site().resolve_region(source, structure_index=2)


def test_policy_context_is_preserved_and_has_no_coordinate_queries(source, monkeypatch):
    original = msm.get

    def guarded(*args, **kwargs):
        assert not kwargs.get('coordinates')
        return original(*args, **kwargs)

    monkeypatch.setattr(msm, 'get', guarded)
    results = []
    for units in (['nm', 'ps'], ['pm', 'fs']):
        with puw.context(standard_units=units):
            results.append(
                site().resolve_region(source, structure_index=2)['projection']
            )
    # The distinct declaration records have the same requested axes and context.
    assert results[0] == results[1]


@pytest.mark.parametrize('form', ['pint', 'unyt'])
def test_application_unit_form_and_policy_are_preserved(source, form):
    with puw.context(default_form=form, standard_units=['pm', 'fs']):
        before = puw.configure.report()
        result = site().resolve_region(source, structure_index=2)
        assert bmt.BindingSite.region_from_dict(result, experimental=True) == result
        assert puw.configure.report() == before


def identity_map(n=6):
    return dict(
        schema='molsysmt.atom_map@1',
        n_source_atoms=n,
        n_output_atoms=n,
        source_to_output_atom_indices=list(range(n)),
        output_to_source_atom_indices=list(range(n)),
    )


def test_unwitnessed_empty_group_requires_explicit_correspondence(source):
    source.topology.atoms['group_index'] = pd.array([0, 1, 0, 1, 1, 1], dtype='Int64')
    target = source.copy()
    options = dict(
        target=target,
        atom_map=identity_map(),
        structure_correspondence=[0, 1, 2],
        structure_index=2,
        expansion='source_members',
    )
    declaration = site('group', [2])
    with pytest.raises(ArgumentError):
        declaration.resolve_region(source, **options, coverage_policy='allow_partial')
    kept = declaration.resolve_region(source, **options, group_correspondence=[0, 1, 2])
    omitted = declaration.resolve_region(
        source,
        **options,
        group_correspondence=[0, 1, -1],
        coverage_policy='allow_partial',
    )
    assert kept['projection']['entries'][0]['state'] == 'empty'
    assert kept['projection']['coverage'] == 'complete'
    assert omitted['projection']['entries'][0]['state'] == 'omitted'
    assert omitted['projection']['coverage'] == 'partial'
    assert bmt.BindingSite.region_from_dict(omitted, experimental=True) == omitted


@pytest.mark.parametrize(
    'vector', [[0, 1, 2], [1, 0, 0], [1, 0], [1, 0, 3], [1, 0, -2], [1, 0, True]]
)
def test_wrong_explicit_group_correspondence_refuses(source, vector):
    with pytest.raises(ArgumentError):
        site('group', [0]).resolve_region(
            source,
            **derived(source),
            group_correspondence=vector,
            expansion='source_members',
            coverage_policy='allow_partial',
        )


@pytest.mark.parametrize('parents', [[0, 1, 2, 2, 1, 2], [0, 0, 0, 2, 0, 2]])
def test_split_or_merged_groups_are_an_unsupported_profile(source, parents):
    target = source.copy()
    target.topology.atoms['group_index'] = pd.array(parents, dtype='Int64')
    with pytest.raises(CapabilityMismatchError):
        site('group', [0]).resolve_region(
            source,
            target=target,
            atom_map=identity_map(),
            structure_correspondence=[0, 1, 2],
            structure_index=2,
            expansion='source_members',
        )


def test_complete_omission_refuses_before_any_indexed_target_query(source, monkeypatch):
    options = derived(source)
    original = msm.get
    indexed_target_queries = []

    def guarded(system, **kwargs):
        if system is options['target'] and 'selection' in kwargs:
            indexed_target_queries.append(kwargs['selection'])
        return original(system, **kwargs)

    monkeypatch.setattr(msm, 'get', guarded)
    declaration = site('group', [0])
    with pytest.raises(ArgumentError):
        declaration.resolve_region(source, **options, expansion='target_group_members')
    assert indexed_target_queries == []


def test_all_omitted_atom_entries_still_keep_the_original_request(source):
    result = site(indices=[0, 3, 4]).resolve_region(
        source, **derived(source), coverage_policy='allow_partial'
    )
    assert result['projection']['exported_target_atom_indices'] == []
    assert result['projection']['coverage'] == 'partial'
    assert [entry['source_index'] for entry in result['projection']['entries']] == [
        0,
        3,
        4,
    ]
    assert result['site']['definition']['binding_support'] == 'unknown'


@pytest.mark.parametrize('vector', [[2, 0], [2, 0, 3], [2, 0, -1], [2, 0, True]])
def test_invalid_complete_structure_relation_refuses(source, vector):
    options = derived(source)
    options['structure_correspondence'] = vector
    with pytest.raises(ArgumentError):
        site().resolve_region(source, **options)


def test_known_declaration_state_and_selected_target_state_must_agree(source):
    with pytest.raises(ArgumentError):
        site(structure_indices=[2], chemical_state_indices=[1]).resolve_region(source)
    options = derived(source)
    msm.set(options['target'], structure_indices=[0], structure_chemical_state_index=1)
    with pytest.raises(ArgumentError):
        site().resolve_region(source, **options)


def test_unsupported_annotation_and_geometry_views_do_not_access_source(monkeypatch):
    declaration = bmt.BindingSite.from_definition(
        source=SOURCE, claim='Hypothesis', criterion={'kind': 'authored'}
    )

    def forbidden(*args, **kwargs):
        pytest.fail('Unsupported view accessed a molecular source')

    monkeypatch.setattr(msm, 'get', forbidden)
    monkeypatch.setattr(msm, 'get_form', forbidden)
    with pytest.raises(CapabilityMismatchError):
        declaration.resolve_region(object())


def test_too_large_map_refuses_before_owner_allocation(source, monkeypatch):
    options = derived(source)
    options['atom_map'] = dict(
        schema='molsysmt.atom_map@1', n_source_atoms=100001, n_output_atoms=3
    )

    def forbidden(*args, **kwargs):
        pytest.fail('Oversized map reached owner normalization')

    monkeypatch.setattr(msm, 'compose_atom_maps', forbidden)
    with pytest.raises(ArgumentError):
        site().resolve_region(source, **options)


def test_invalid_owner_atom_map_is_not_accepted_by_resealing(source):
    record = site().resolve_region(source, **derived(source))
    record['correspondence']['atom_map']['source_to_output_atom_indices'][5] = 1
    record['content_sha256'] = _digest(record)
    with pytest.raises(StructuralInconsistencyError):
        bmt.BindingSite.region_from_dict(record, experimental=True)


@pytest.mark.parametrize(
    'field,value',
    [
        ('format_version', '1.0'),
        ('experimental', False),
        ('record_id', 'wrong'),
        ('content_sha256', '0' * 64),
        ('producer', {}),
        ('unexpected', 'not an admitted extension'),
    ],
)
def test_result_reader_is_strict_about_version_fields_and_integrity(
    source, field, value
):
    record = site().resolve_region(source, structure_index=2)
    record[field] = deepcopy(value)
    if field != 'content_sha256':
        record['content_sha256'] = _digest(record)
    with pytest.raises((ArgumentError, CapabilityMismatchError)):
        bmt.BindingSite.region_from_dict(record, experimental=True)


@pytest.mark.parametrize('field', ['atom_map', 'canonical_atom_map'])
def test_missing_required_map_material_refuses_in_historical_target_record(
    source, field
):
    record = site().resolve_region(source, **derived(source))
    record['correspondence'][field] = None
    record['content_sha256'] = _digest(record)
    with pytest.raises(ArgumentError):
        bmt.BindingSite.region_from_dict(record, experimental=True)


def test_resource_failure_during_data_only_owner_admission_propagates(
    source, monkeypatch
):
    record = site().resolve_region(source, structure_index=2)
    failure = OSError(errno.ENOMEM, 'owner map allocation')

    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(msm, 'compose_atom_maps', fail)
    with pytest.raises(OSError) as caught:
        bmt.BindingSite.region_from_dict(record, experimental=True)
    assert caught.value is failure
