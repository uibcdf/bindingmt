"""Malformed saved operation options keep the public diagnostic boundary."""

import json
from copy import deepcopy
from pathlib import Path

import molsysmt as msm
import pytest
import pyunitwizard as puw
from pyunitwizard.record import QuantityRecord

import bindingmt as bmt
from bindingmt._private.interface_records import _digest


@pytest.fixture(scope='module', params=[1, 2], ids=['positive', 'empty'])
def records(request):
    import bmt_test_support

    data = Path(bmt_test_support.__path__[0]) / 'site_interface/source.h5msm'
    if not data.is_file():
        data = (
            Path(__file__).resolve().parents[1] / 'examples/site_interface/source.h5msm'
        )
    before = data.read_bytes()
    source = msm.h5msm.read(str(data))
    result = bmt.calculate_interface_contacts(
        source,
        participants={'region': [2], 'partner': [3]},
        cutoff='0.3 nm',
        structure_index=request.param,
        source_id='authored-options-control',
        max_dense_pairs=16,
    )
    assert result.contact_calculation['local_pairs'] == (
        [[2, 3]] if request.param == 1 else []
    )
    footprint = result.project_footprint(participant='region')
    yield {'interface': result.to_dict(), 'footprint': footprint.to_dict()}
    assert data.read_bytes() == before


def block_molecular(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Saved admission attempted molecular work')

    for name in ('select', 'get', 'set', 'convert', 'extract', 'copy'):
        monkeypatch.setattr(msm, name, forbidden)
    monkeypatch.setattr(msm.h5msm, 'read', forbidden)
    monkeypatch.setattr(msm.structure, 'get_contacts', forbidden)
    monkeypatch.setattr(msm.structure, 'get_distances', forbidden)


def read_record(kind, record, route, tmp_path):
    cls = bmt.BindingInterface if kind == 'interface' else bmt.ContactFootprint
    if route == 'dict':
        return cls.from_dict(record, experimental=True)
    return cls.read(tmp_path / 'record.json', experimental=True)


@pytest.mark.parametrize('kind', ('interface', 'footprint'))
@pytest.mark.parametrize('position', (0, 1))
@pytest.mark.parametrize(
    'invalid', (None, False, 17, 'invalid', [], [{}], {}, 'missing')
)
@pytest.mark.parametrize('route', ('dict', 'file'))
def test_malformed_contact_options_use_calls_diagnostic(
    records, kind, position, invalid, route, tmp_path, monkeypatch
):
    record = deepcopy(records[kind])
    calculation = (
        record['analyses']['contacts']['calculation']
        if kind == 'interface'
        else record['calculation']
    )
    if invalid == 'missing':
        del calculation['calls'][position]['options']
    else:
        calculation['calls'][position]['options'] = deepcopy(invalid)
    # A valid outer seal must not mask the actual structural admission check.
    record['content_sha256'] = _digest(record)
    before = deepcopy(record)
    path = tmp_path / 'record.json'
    payload = json.dumps(record).encode()
    path.write_bytes(payload)
    block_molecular(monkeypatch)
    with pytest.raises(bmt.exceptions.ArgumentError) as caught:
        read_record(kind, record, route, tmp_path)
    expected_field = 'participants' if position == 0 and invalid == {} else 'calls'
    assert caught.value.extra['arg_name'] == expected_field
    assert record == before
    assert path.read_bytes() == payload


@pytest.mark.parametrize('kind', ('interface', 'footprint'))
@pytest.mark.parametrize(
    'form,units',
    [('pint', ['nm', 'ps']), ('pint', ['pm', 'fs']), ('unyt', ['pm', 'fs'])],
)
def test_valid_saved_contact_options_keep_original_evidence(
    records, kind, form, units, tmp_path, monkeypatch
):
    record = records[kind]
    payload = json.dumps(record).encode()
    (tmp_path / 'record.json').write_bytes(payload)
    block_molecular(monkeypatch)
    with puw.context(default_form=form, default_parser='pint', standard_units=units):
        policy = puw.configure.report()
        for route in ('dict', 'file'):
            restored = read_record(kind, record, route, tmp_path)
            assert restored.to_dict() == record
            assert restored.copy().to_dict() == record
            summary = restored.summary()
            if kind == 'footprint':
                assert summary['current_assessment'] == 'not_performed'
            else:
                assert (
                    summary['analyses']['contacts']['original_outcome'] == 'completed'
                )
        assert puw.configure.report() == policy
    assert (tmp_path / 'record.json').read_bytes() == payload


@pytest.mark.parametrize('kind', ('interface', 'footprint'))
@pytest.mark.parametrize('route', ('dict', 'file'))
@pytest.mark.parametrize(
    'error',
    (
        MemoryError('controlled'),
        OSError(12, 'controlled'),
        KeyboardInterrupt('controlled'),
        SystemExit(7),
    ),
)
def test_quantity_provider_resource_and_process_errors_keep_identity(
    records, kind, route, error, tmp_path, monkeypatch
):
    record = records[kind]
    payload = json.dumps(record).encode()
    (tmp_path / 'record.json').write_bytes(payload)

    def fail(*args, **kwargs):
        raise error

    block_molecular(monkeypatch)
    monkeypatch.setattr(QuantityRecord, 'from_dict', fail)
    with pytest.raises(type(error)) as caught:
        read_record(kind, record, route, tmp_path)
    assert caught.value is error
    assert (tmp_path / 'record.json').read_bytes() == payload
