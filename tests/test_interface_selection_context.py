"""Participant queries use the actual chosen frame/state, never defaults."""

import errno
import hashlib
from copy import deepcopy
from pathlib import Path

import molsysmt as msm
import numpy as np
import pytest
import pyunitwizard as puw
from molsysmt._private.smonitor import StructuralInconsistencyError
from molsysmt.native import Structures

import bindingmt as bmt

INPUT = Path(__file__).parent / 'data/interface_native_handoff/cb6af716/source.h5msm'
QUERIES = {'A': 'formal_charge>0', 'B': 'formal_charge<0'}


@pytest.fixture
def system():
    """An authored derivative; the historical four-ion asset is never rewritten."""
    original = INPUT.read_bytes()
    source = msm.h5msm.read(str(INPUT))
    source.interactions = {}
    source.chemical_states.append_state()
    states = msm.convert(source.chemical_states, to_form='molsysmt.ChemicalStatesDict')
    # The newly authored monatomic state has complete, empty connectivity.
    states.data['states'][1]['connectivity_completeness'] = 'complete'
    source.chemical_states = msm.convert(states, to_form='molsysmt.ChemicalStates')
    msm.set(source, chemical_state=0, formal_charge=[-1, 1, 0, 0])
    msm.set(source, chemical_state=1, formal_charge=[0, 0, 1, -1])
    points = [
        [[0.25, 0, 0], [0, 0, 0], [5, 0, 0], [9, 0, 0]],
        [[5, 0, 0], [9, 0, 0], [0, 0, 0], [0.25, 0, 0]],
    ]
    source.structures = Structures(coordinates=puw.quantity(points, 'nm'))
    msm.set(source, structure_chemical_state_index=[0, 1])
    source.interactions = {
        f'ionic-{frame}': msm.interactions.ionic.get_ionic_interactions(
            source,
            '0.3 nm',
            selection='all',
            structure_indices=[frame],
            chemical_state='structure',
            selection_mode='internal',
            pbc=False,
        )
        for frame in (0, 1)
    }
    yield source
    assert INPUT.read_bytes() == original


def calculate(system, method, *, participants=None, frame=1):
    options = {
        'participants': QUERIES if participants is None else participants,
        'structure_index': frame,
        'source_id': 'authored-state-switch',
    }
    if method == 'contacts':
        return bmt.calculate_interface_contacts(system, cutoff='0.3 nm', **options)
    if method == 'burial':
        return bmt.calculate_interface_burial(
            system,
            probe_radius='0.14 nm',
            n_sphere_points=120,
            hydrogen_policy='as_supplied',
            **options,
        )
    return bmt.characterize_interface_interactions(
        system,
        grouping=[['A'], ['B']],
        analysis_names=[f'ionic-{frame}'],
        **options,
    )


def membership(result):
    return {k: v['source_atom_indices'] for k, v in result.participants.items()}


@pytest.mark.parametrize('method', ['contacts', 'burial', 'ionic'])
@pytest.mark.parametrize('query', ['chemical', 'spatial', 'combined'])
def test_selected_frame_and_state_resolve_actual_participants(
    system, monkeypatch, method, query
):
    expected = {'A': [2], 'B': [3]}
    participants = QUERIES
    if query == 'spatial':
        participants = {'A': 'all within 0.3 nm of atom_index==2', 'B': [1]}
        expected = {'A': [2, 3], 'B': [1]}
    elif query == 'combined':
        participants = {
            'A': 'formal_charge>0 within 0.3 nm of atom_index==3',
            'B': 'formal_charge<0',
        }
    # Independent fixed expectations, and the actual owner query as a control.
    for name, selection in participants.items():
        assert (
            msm.select(
                system,
                selection=selection,
                structure_indices=[1],
                chemical_state='structure',
            )
            == expected[name]
        )
    assert msm.select(system, selection='formal_charge>0') == [1]
    assert msm.select(system, selection='formal_charge<0') == [0]
    positions = puw.get_value(system.structures.coordinates, to_unit='nm').copy()
    charges = [
        list(msm.get(system, chemical_state=state, formal_charge=True))
        for state in (0, 1)
    ]
    calls = []
    original = msm.select

    def traced(source, selection='all', **options):
        if (
            source is system
            and isinstance(selection, str)
            and selection
            in (value for value in participants.values() if isinstance(value, str))
        ):
            calls.append((selection, deepcopy(options)))
        return original(source, selection=selection, **options)

    monkeypatch.setattr(msm, 'select', traced)
    result = calculate(system, method, participants=participants)
    assert membership(result) == expected
    assert len(calls) == sum(isinstance(v, str) for v in participants.values())
    for _, options in calls:
        assert options['structure_indices'] == [1]
        assert options['chemical_state'] == 'structure'
    if method == 'contacts':
        assert result.contact_calculation['local_pairs'] == (
            [] if query == 'spatial' else [[2, 3]]
        )
    elif method == 'burial':
        areas = [
            float(puw.get_value(a, to_unit='nm**2'))
            for a in result.burial_areas['sides'].values()
        ]
        assert all(a == 0 if query == 'spatial' else a > 0 for a in areas)
        assert (
            result.surface_calculation['context']['structure_chemical_state_index'] == 1
        )
    else:
        view = result.interaction_analyses['ionic-1']
        assert len(view['observations']) == (0 if query == 'spatial' else 1)
        assert view['context']['structure_index'] == 1
        assert view['context']['chemical_state_index'] == 1
    np.testing.assert_array_equal(
        puw.get_value(system.structures.coordinates, to_unit='nm'), positions
    )
    assert charges == [
        list(msm.get(system, chemical_state=state, formal_charge=True))
        for state in (0, 1)
    ]


@pytest.mark.parametrize('method', ['contacts', 'burial', 'ionic'])
@pytest.mark.parametrize('frame', [0, 1])
def test_explicit_numeric_members_keep_their_chosen_frame(system, method, frame):
    expected = {'A': [1], 'B': [0]} if frame == 0 else {'A': [2], 'B': [3]}
    assert (
        membership(calculate(system, method, participants=expected, frame=frame))
        == expected
    )


@pytest.mark.parametrize('method', ['contacts', 'burial', 'ionic'])
def test_unknown_association_refuses_before_required_calculation(
    system, monkeypatch, method
):
    msm.set(system, structure_indices=[1], structure_chemical_state_index=[None])

    def forbidden(*args, **kwargs):
        pytest.fail('Unresolved selection reached a numerical/scoped-observation child')

    monkeypatch.setattr(msm.structure, 'get_contacts', forbidden)
    monkeypatch.setattr(msm.structure, 'get_distances', forbidden)
    monkeypatch.setattr(msm.physchem, 'get_sasa', forbidden)
    monkeypatch.setattr(msm.Interactions, 'between_selections', forbidden)
    with pytest.raises(StructuralInconsistencyError):
        calculate(system, method)


@pytest.mark.parametrize('method', ['contacts', 'burial', 'ionic'])
@pytest.mark.parametrize('kind', ['memory', 'enomem', 'interrupt', 'ordinary'])
def test_selected_participant_query_preserves_owner_failure_identity(
    system, monkeypatch, method, kind
):
    error = {
        'memory': MemoryError('owned selection allocation'),
        'enomem': OSError(errno.ENOMEM, 'owned selection allocation'),
        'interrupt': KeyboardInterrupt('owned selection stop'),
        'ordinary': ValueError('owned selection rejection'),
    }[kind]

    def failed(source, selection='all', **options):
        assert source is system
        assert options['structure_indices'] == [1]
        assert options['chemical_state'] == 'structure'
        raise error

    monkeypatch.setattr(msm, 'select', failed)
    with pytest.raises(type(error)) as caught:
        calculate(system, method)
    assert caught.value is error


@pytest.mark.parametrize('method', ['contacts', 'burial', 'ionic'])
def test_selected_result_history_reads_without_molecular_reselection(
    system, monkeypatch, method, tmp_path
):
    result = calculate(system, method)
    assert membership(result) == {'A': [2], 'B': [3]}
    path = tmp_path / 'selected.json'
    result.save(path)
    before = path.read_bytes()
    digest = hashlib.sha256(before).hexdigest()

    def forbidden(*args, **kwargs):
        pytest.fail('Historical result attempted molecular reassessment')

    for name in ['select', 'get', 'convert', 'extract']:
        monkeypatch.setattr(msm, name, forbidden)
    with puw.context(default_form='unyt', default_parser='pint'):
        restored = bmt.BindingInterface.read(path, experimental=True)
        assert restored.to_dict() == result.to_dict()
    assert path.read_bytes() == before
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_selected_ionic_context_composes_with_actual_contact_surface_results(system):
    original = calculate(system, 'ionic')
    children = {name: calculate(system, name) for name in ['contacts', 'burial']}
    annotated = original.with_characterizations(
        system, characterizations=children, participant='A'
    )
    assert annotated.analyses == original.analyses
    for name, child in children.items():
        assert annotated.characterizations[name]['record'] == child.to_dict()
        assert membership(child) == {'A': [2], 'B': [3]}
