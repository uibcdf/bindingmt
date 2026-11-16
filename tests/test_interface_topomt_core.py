"""Frozen actual owner-derived records remain readable without optional owners."""

import builtins
import json
import socket
from copy import deepcopy

import pytest
import pyunitwizard as puw
from pyunitwizard.record import RecordError

import bindingmt as bmt
from bindingmt._private.interface_records import _digest
from bmt_test_support.pisa_case import RECORDS, verify_records

FILES = [
    'I-positive.json',
    'I-empty.json',
    'I-failed.json',
    'I-unknown.json',
    'F-A.json',
    'F-B.json',
    'F-empty.json',
    'F-mixed.json',
]


def test_absent_selected_owner_has_truthful_installation_advice(monkeypatch):
    from depdigest.core import checker

    installed = checker.is_installed
    monkeypatch.setattr(
        checker,
        'is_installed',
        lambda name: False if name == 'topomt' else installed(name),
    )
    with pytest.raises(bmt.exceptions.LibraryNotFoundError) as caught:
        bmt.BindingInterface.from_topomt(
            None,
            question={'kind': 'association'},
            participants={'A': {'id': 'a'}, 'B': {'id': 'b'}},
            grouping=[['A'], ['B']],
        )
    assert 'conda install -c uibcdf topomt' in str(caught.value)
    assert 'pip install topomt' not in str(caught.value)


def test_selected_owner_transitive_import_error_is_preserved(monkeypatch):
    from depdigest.core import checker

    installed = checker.is_installed
    monkeypatch.setattr(
        checker,
        'is_installed',
        lambda name: True if name == 'topomt' else installed(name),
    )
    original = builtins.__import__
    error = ModuleNotFoundError('Controlled missing transitive dependency')

    def fail(module, *args, **kwargs):
        if module == 'topomt.provider_output':
            raise error
        return original(module, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', fail)
    with pytest.raises(ModuleNotFoundError) as caught:
        bmt.BindingInterface.from_topomt(
            None,
            question={'kind': 'association'},
            participants={'A': {'id': 'a'}, 'B': {'id': 'b'}},
            grouping=[['A'], ['B']],
        )
    assert caught.value is error


def read(name, path=None):
    cls = bmt.BindingInterface if name.startswith('I-') else bmt.ContactFootprint
    return cls.read(RECORDS / name if path is None else path, experimental=True)


def reverse_objects(value):
    if isinstance(value, dict):
        return {
            key: reverse_objects(item) for key, item in reversed(list(value.items()))
        }
    if isinstance(value, list):
        return [reverse_objects(item) for item in value]
    return value


@pytest.mark.parametrize('name', FILES)
@pytest.mark.parametrize(
    'form,units',
    [('pint', ['nm', 'ps']), ('pint', ['pm', 'fs']), ('unyt', ['pm', 'fs'])],
)
def test_frozen_original_records_read_with_no_optional_import_or_network(
    name, form, units, monkeypatch, tmp_path
):
    manifest = verify_records()
    payload = (RECORDS / name).read_bytes()
    original = builtins.__import__

    def guarded(module, *args, **kwargs):
        if module.split('.')[0] in {
            'topomt',
            'pharmacophoremt',
            'pisa',
            'vina',
            'molsysviewer',
        }:
            raise AssertionError('Core read imported an optional owner/engine')
        return original(module, *args, **kwargs)

    def forbidden(*args, **kwargs):
        raise AssertionError('Core read attempted network access')

    monkeypatch.setattr(builtins, '__import__', guarded)
    monkeypatch.setattr(socket.socket, 'connect', forbidden)
    with puw.context(default_form=form, default_parser='pint', standard_units=units):
        before = puw.configure.report()
        result = read(name)
        record = json.loads(payload)
        assert result.to_dict() == record == result.copy().to_dict()
        assert result.summary()['current_assessment'] == 'not_performed'
        result.save(tmp_path / name)
        assert read(name, tmp_path / name).to_dict() == record
        assert puw.configure.report() == before
        if name.startswith('F-'):
            expected = (
                []
                if name == 'F-empty.json'
                else manifest['independent_expectations']['other_bonds'][
                    'B' if name == 'F-B.json' else 'A'
                ]
            )
            assert result.source_atom_indices == result.local_atom_indices == expected
        else:
            outcome = (
                'completed'
                if name in ('I-positive.json', 'I-empty.json')
                else 'failed'
                if name == 'I-failed.json'
                else 'unknown'
            )
            assert result.summary()['analyses']['pisa']['original_outcome'] == outcome
    assert (RECORDS / name).read_bytes() == payload


@pytest.mark.parametrize('name', FILES)
def test_reversed_or_sorted_json_keeps_meaningful_order_and_original_digests(
    name, tmp_path
):
    original = read(name).to_dict()
    for index, altered in enumerate(
        [reverse_objects(original), json.loads(json.dumps(original, sort_keys=True))]
    ):
        path = tmp_path / f'{index}.json'
        path.write_text(json.dumps(altered))
        restored = read(name, path)
        assert restored.to_dict() == original
        if name == 'F-mixed.json':
            assert [g['name'] for g in restored.summary()['contact_groups']] == [
                'h-bonds',
                'other-bonds',
            ]


@pytest.mark.parametrize(
    'mutation',
    [
        'parent',
        'member',
        'witness',
        'local_axis',
        'state_unknown',
        'state_changed',
        'geometry',
        'scope',
        'current_receipt',
        'distances',
        'group_order',
        'receiver',
        'bool_index',
        'bool_map',
        'bool_witness',
        'float_index',
        'counterpart_chain',
        'submitted_axis',
    ],
)
def test_resealed_inconsistent_reported_footprints_refuse(mutation):
    record = read('F-A.json').to_dict()
    projection = record['reported_projection']
    if mutation == 'parent':
        projection['parent_ref']['record_id'] = 'different'
    elif mutation == 'member':
        record['footprint']['source_atom_indices'][0] = 4
    elif mutation == 'witness':
        record['footprint']['contributors'][0] = []
    elif mutation == 'local_axis':
        projection['local_atom_map'][0] = [0, 1]
    elif mutation == 'state_unknown':
        projection['observation']['metadata']['chemical_state_context'] = None
    elif mutation == 'state_changed':
        projection['source_context']['chemical_state_index'] += 1
    elif mutation == 'geometry':
        projection['source_context']['geometry_sha256_nm_float64'] = None
    elif mutation == 'scope':
        projection['source_context']['scope'] = 'current_finding'
    elif mutation == 'current_receipt':
        projection['source_context']['assessment_kind'] = 'current_guarantee'
    elif mutation == 'distances':
        projection['reported_distances']['other-bonds'] = projection['observation'][
            'findings'
        ][0]['measurements']['interface_area']['value']
    elif mutation == 'group_order':
        projection['contact_groups'] = []
    elif mutation == 'receiver':
        projection['receiver']['original_execution'] = True
    elif mutation == 'bool_map':
        projection['local_atom_map'][0] = [False, False]
    elif mutation == 'bool_witness':
        record['footprint']['contributors'][0][0]['endpoint_index'] = False
    elif mutation == 'float_index':
        record['footprint']['source_atom_indices'][0] = 0.0
    elif mutation == 'counterpart_chain':
        projection['observation']['findings'][0]['fields']['contacts']['other-bonds'][
            0
        ]['original']['chain-2'] = 'A'
    elif mutation == 'submitted_axis':
        projection['observation']['metadata']['selected_atom_indices'].append(40)
    else:
        record['footprint']['source_atom_indices'][0] = False
    record['content_sha256'] = _digest(record)
    before = deepcopy(record)
    with pytest.raises((bmt.exceptions.ArgumentError, RecordError)):
        bmt.ContactFootprint.from_dict(record, experimental=True)
    assert record == before


@pytest.mark.parametrize('name', FILES)
def test_frozen_records_need_opt_in_and_exclusive_output(name, tmp_path):
    cls = bmt.BindingInterface if name.startswith('I-') else bmt.ContactFootprint
    with pytest.raises(bmt.exceptions.ArgumentError):
        cls.read(RECORDS / name)
    result = read(name)
    path = tmp_path / name
    result.save(path)
    payload = path.read_bytes()
    with pytest.raises(FileExistsError):
        result.save(path)
    assert path.read_bytes() == payload
