"""Actual optional owner intake and original-source reported footprint receiving."""

import hashlib
import importlib.util
import socket
import subprocess
import zipfile
from copy import deepcopy

import molsysmt as msm
import pytest
import pyunitwizard as puw
from argdigest import UnknownArgumentError

import bindingmt as bmt
from bindingmt._private.interface_records import _digest
from bindingmt._private.interface_reported_contacts import aggregate_endpoints
from bindingmt._private.smonitor import CapabilityMismatchError
from bmt_test_support.pisa_case import ORIGINAL, source, verify


@pytest.fixture(scope='module')
def tmt():
    if importlib.util.find_spec('topomt') is None:
        pytest.skip('The optional TopoMT saved-result owner is not installed.')
    import topomt

    return topomt


@pytest.fixture(scope='module')
def cases(tmt):
    from topomt.provider_output import ProviderRun

    positive = source()
    empty = source(empty=True)
    return {
        'positive': tmt.load_provider_output(
            ORIGINAL / 'pisa.zip', molecular_system=positive
        ),
        'empty': tmt.load_provider_output(
            ORIGINAL / 'pisa-empty.zip', molecular_system=empty
        ),
        'raw': ProviderRun.load(ORIGINAL / 'pisa.zip'),
        'failed': ProviderRun.load(ORIGINAL / 'pisa-failed.zip'),
        'source': positive,
        'empty_source': empty,
    }


def intake(value, **kwargs):
    return bmt.BindingInterface.from_topomt(
        value,
        question={'kind': 'association'},
        participants={
            'A': {
                'source_id': 'original-pisa-control',
                'source_atom_indices': list(range(20)),
            },
            'B': {
                'source_id': 'original-pisa-control',
                'source_atom_indices': list(range(20, 40)),
            },
        },
        grouping=[['A'], ['B']],
        **kwargs,
    )


@pytest.fixture(scope='module')
def interface(cases):
    return intake(cases['positive'])


def project(
    interface,
    molecular_system,
    *,
    participant='A',
    groups=('other-bonds',),
    bundle_file=None,
    **kwargs,
):
    return interface.project_reported_footprint(
        molecular_system,
        bundle_file=ORIGINAL / 'pisa.zip' if bundle_file is None else bundle_file,
        finding_id='1',
        participant=participant,
        contact_groups=groups,
        **kwargs,
    )


def forbidden(*args, **kwargs):
    raise AssertionError('Unexpected molecular/original-engine/process/network work')


def block_molecular(monkeypatch):
    for name in ('get', 'set', 'select', 'copy', 'convert', 'extract'):
        monkeypatch.setattr(msm, name, forbidden)
    monkeypatch.setattr(msm.h5msm, 'read', forbidden)
    monkeypatch.setattr(msm.structure, 'get_contacts', forbidden)
    monkeypatch.setattr(msm.physchem, 'get_sasa', forbidden)


@pytest.mark.parametrize(
    'case,outcome,count',
    [
        ('positive', 'completed', 1),
        ('empty', 'completed', 0),
        ('failed', 'failed', None),
        ('raw', 'unknown', None),
    ],
)
def test_actual_original_outcomes_and_independent_roundtrip(
    cases, tmt, case, outcome, count, tmp_path, monkeypatch
):
    original = cases[case]
    run = original.run if hasattr(original, 'run') else original
    before = (run.metadata, run.attribution, run.artifacts)
    block_molecular(monkeypatch)
    monkeypatch.setattr(tmt, 'load_provider_output', forbidden)
    monkeypatch.setattr(subprocess, 'Popen', forbidden)
    monkeypatch.setattr(socket.socket, 'connect', forbidden)
    result = intake(
        original, artifact_ref={'uri': 'unavailable-original.zip', 'sha256': 'a' * 64}
    )
    summary = result.summary()
    assert summary['analyses']['pisa']['original_outcome'] == outcome
    assert summary['analyses']['pisa']['finding_count'] == count
    view = result.topomt_analyses['pisa']
    assert view['metadata'] == run.metadata and view['attribution'] == run.attribution
    assert view['artifact_checksums'] == {
        name: hashlib.sha256(data).hexdigest() for name, data in run.artifacts
    }
    assert (run.metadata, run.attribution, run.artifacts) == before
    path = tmp_path / 'I.json'
    result.save(path)
    restored = bmt.BindingInterface.read(path, experimental=True)
    assert restored.to_dict() == result.to_dict() == restored.copy().to_dict()
    assert restored.summary()['current_assessment'] == 'not_performed'
    if outcome == 'failed':
        assert restored.analyses['pisa']['failure']['stage'] == 'analysis'
    receiver = restored.analyses['pisa']['extensions']['bindingmt.original_intake']
    assert receiver['original_execution'] is False


def test_unsupplied_analysis_is_distinct_from_completed_zero(cases):
    absent = bmt.characterize_interface(
        question={'kind': 'association'},
        participants={'A': {'id': 'a'}, 'B': {'id': 'b'}},
        grouping=[['A'], ['B']],
    )
    assert absent.summary()['analysis_presence'] == 'not_supplied'
    assert intake(cases['empty']).summary()['analysis_presence'] == 'supplied'


@pytest.mark.parametrize(
    'form,units',
    [('pint', ['nm', 'ps']), ('pint', ['pm', 'fs']), ('unyt', ['pm', 'fs'])],
)
def test_original_quantities_keep_explicit_units_and_application_policy(
    cases, form, units, tmp_path
):
    with puw.context(default_form=form, default_parser='pint', standard_units=units):
        before = puw.configure.report()
        result = intake(cases['positive'])
        original = result.topomt_analyses['pisa']['findings'][0]
        from pyunitwizard.record import QuantityRecord

        value = QuantityRecord.from_dict(
            original['measurements']['interface_area']['value']
        ).to_quantity(
            field='bindingmt.topomt.pisa.interface.1.interface_area',
            unit='nm**2',
            dimensionality={'[L]': 2},
        )
        assert float(puw.get_value(value, to_unit='nm**2')) == pytest.approx(
            1.2573203716, rel=2e-14
        )
        result.save(tmp_path / 'I.json')
        assert (
            bmt.BindingInterface.read(tmp_path / 'I.json', experimental=True).to_dict()
            == result.to_dict()
        )
        assert puw.configure.report() == before


@pytest.mark.parametrize('participant', ['A', 'B'])
@pytest.mark.parametrize(
    'groups',
    [
        ('other-bonds',),
        ('h-bonds',),
        ('h-bonds', 'other-bonds'),
        ('other-bonds', 'h-bonds'),
    ],
)
def test_actual_reported_members_contributors_and_saved_footprint(
    interface, cases, participant, groups, tmp_path
):
    before = interface.to_dict()
    result = project(interface, cases['source'], participant=participant, groups=groups)
    expected = (
        verify()['expected']['other_bonds'][participant]
        if 'other-bonds' in groups
        else []
    )
    assert result.source_atom_indices == result.local_atom_indices == expected
    assert result.summary()['value_state'] == (
        'evaluated_nonempty' if expected else 'evaluated_empty'
    )
    projection = result.to_dict()['reported_projection']
    assert projection['contact_groups'] == list(groups)
    assert projection['parent_ref'] == {
        'record_id': interface.record_id,
        'content_sha256': before['content_sha256'],
    }
    assert projection['source_context']['assessment_kind'] == 'historical_check_receipt'
    for member, witnesses in zip(
        result.source_atom_indices, result.to_dict()['footprint']['contributors']
    ):
        side = 0 if participant == 'A' else 1
        independently_expected = [
            {'group': group, 'contact_index': i, 'endpoint_index': side}
            for group in groups
            for i, row in enumerate(
                projection['observation']['findings'][0]['fields']['contacts'][group]
            )
            if row['source_atom_indices'][side] == member
        ]
        assert witnesses == independently_expected
    result.save(tmp_path / 'F.json')
    restored = bmt.ContactFootprint.read(tmp_path / 'F.json', experimental=True)
    assert restored.to_dict() == result.to_dict() == result.copy().to_dict()
    assert (
        restored.check_context(cases['source'], bundle_file=ORIGINAL / 'pisa.zip')[
            'result'
        ]
        == 'match'
    )
    assert restored.summary()['current_assessment'] == 'not_performed'
    assert interface.to_dict() == before


def test_nonempty_overlapping_groups_keep_canonical_members_aligned_map_and_every_witness():
    # Independent test-local group IDs and explicit original occurrence indices.
    g1 = [
        {'source_atom_index': 7, 'contact_index': 0, 'endpoint_index': 0},
        {'source_atom_index': 2, 'contact_index': 1, 'endpoint_index': 0},
        {'source_atom_index': 7, 'contact_index': 2, 'endpoint_index': 0},
    ]
    g2 = [
        {'source_atom_index': 2, 'contact_index': 0, 'endpoint_index': 0},
        {'source_atom_index': 11, 'contact_index': 1, 'endpoint_index': 0},
    ]
    expected = [
        [
            {'group': 'g1', 'contact_index': 1, 'endpoint_index': 0},
            {'group': 'g2', 'contact_index': 0, 'endpoint_index': 0},
        ],
        [
            {'group': 'g1', 'contact_index': 0, 'endpoint_index': 0},
            {'group': 'g1', 'contact_index': 2, 'endpoint_index': 0},
        ],
        [{'group': 'g2', 'contact_index': 1, 'endpoint_index': 0}],
    ]
    result = aggregate_endpoints([('g1', g1), ('g2', g2)], {2: 6, 7: 1, 11: 4})
    assert result['source_atom_indices'] == [2, 7, 11]
    assert result['local_atom_indices'] == [6, 1, 4]
    assert result['contributors'] == expected
    reversed_result = aggregate_endpoints(
        [('g2', list(reversed(g2))), ('g1', list(reversed(g1)))], {2: 6, 7: 1, 11: 4}
    )
    assert reversed_result['source_atom_indices'] == [2, 7, 11] and reversed_result[
        'local_atom_indices'
    ] == [6, 1, 4]
    assert [
        sorted(w, key=lambda x: (x['group'], x['contact_index']))
        for w in reversed_result['contributors']
    ] == [sorted(w, key=lambda x: (x['group'], x['contact_index'])) for w in expected]


@pytest.mark.parametrize(
    'groups', [[], ['other-bonds', 'other-bonds'], ['unknown'], 'other-bonds', None]
)
def test_invalid_group_requests_refuse_before_recovery(
    interface, cases, tmt, groups, monkeypatch
):
    monkeypatch.setattr(tmt, 'load_provider_output', forbidden)
    with pytest.raises(bmt.exceptions.ArgumentError):
        project(interface, cases['source'], groups=groups)


@pytest.mark.parametrize('groups', [['salt-bridges'], ['ss-bonds'], ['cov-bonds']])
def test_known_but_unqualified_groups_remain_unsupported(
    interface, cases, tmt, groups, monkeypatch
):
    monkeypatch.setattr(tmt, 'load_provider_output', forbidden)
    with pytest.raises(CapabilityMismatchError):
        project(interface, cases['source'], groups=groups)


@pytest.mark.parametrize(
    'scope', ['current_finding', 'placed_geometry', 'preserved_evidence']
)
def test_stronger_or_other_scope_cannot_silently_downgrade(
    interface, cases, tmt, scope, monkeypatch
):
    monkeypatch.setattr(tmt, 'load_provider_output', forbidden)
    with pytest.raises(CapabilityMismatchError):
        project(interface, cases['source'], assurance_scope=scope)


def test_zero_analysis_unknown_and_original_failure_cannot_become_an_empty_projection(
    cases,
):
    for name in ['empty', 'raw', 'failed']:
        original = intake(cases[name])
        before = original.to_dict()
        with pytest.raises(bmt.exceptions.ArgumentError):
            project(
                original, cases['empty_source'] if name == 'empty' else cases['source']
            )
        assert original.to_dict() == before


@pytest.mark.parametrize('groups', [('other-bonds',), ('h-bonds',)])
def test_unavailable_endpoint_map_never_becomes_evaluated_empty(cases, groups):
    from topomt.third_party.output import ProviderRecord
    from topomt.third_party.pisa.output import PisaOutput

    original = cases['positive'].interfaces[0]
    fields = original.fields
    for rows in fields['contacts'].values():
        for row in rows:
            row.pop('source_atom_indices', None)
    # Explicitly derived software carrier, not a newly recovered original capture.
    unmapped = ProviderRecord(
        original.source_id,
        'interface',
        original.run,
        atom_indices=None,
        atom_role='unavailable_contact_atom_map',
        _fields=fields,
    )
    result = intake(PisaOutput(original.run, (unmapped,)))
    with pytest.raises(bmt.exceptions.ArgumentError) as caught:
        project(result, cases['source'], groups=groups)
    assert caught.value.extra['arg_name'] == 'mapping'


@pytest.mark.parametrize('groups', [('other-bonds',), ('h-bonds',)])
def test_historical_match_survives_read_but_changed_current_source_refuses(
    interface, cases, groups, tmp_path, monkeypatch
):
    footprint = project(interface, cases['source'], groups=groups)
    before = footprint.to_dict()
    footprint.save(tmp_path / 'F.json')
    changed = msm.copy(cases['source'])
    coordinates = msm.get(
        changed, element='atom', selection=[4], structure_indices=[0], coordinates=True
    )
    moved = puw.get_value(coordinates, to_unit='nm').copy()
    moved[0, 0, 0] += 0.01
    msm.set(
        changed,
        selection=[4],
        structure_indices=[0],
        coordinates=puw.quantity(moved, 'nm'),
    )
    with monkeypatch.context() as patch:
        block_molecular(patch)
        restored = bmt.ContactFootprint.read(tmp_path / 'F.json', experimental=True)
        assert (
            restored.to_dict() == before
            and restored.summary()['current_assessment'] == 'not_performed'
        )
    with pytest.raises(ValueError):
        restored.check_context(changed, bundle_file=ORIGINAL / 'pisa.zip')
    with pytest.raises(ValueError):
        project(interface, changed, groups=groups)
    assert restored.to_dict() == before
    assert interface.analyses['pisa']['outcome'] == 'completed'


def test_missing_asset_preserves_positive_or_empty_history(interface, cases, tmp_path):
    for groups in [('other-bonds',), ('h-bonds',)]:
        footprint = project(interface, cases['source'], groups=groups)
        before = footprint.to_dict()
        with pytest.raises(FileNotFoundError):
            footprint.check_context(cases['source'], bundle_file=tmp_path / 'gone.zip')
        assert footprint.to_dict() == before


def test_preserved_exact_asset_digest_is_required_before_owner_recovery(
    cases, tmt, monkeypatch
):
    result = intake(
        cases['positive'], artifact_ref={'uri': 'declared.zip', 'sha256': 'a' * 64}
    )
    before = result.to_dict()
    monkeypatch.setattr(tmt, 'load_provider_output', forbidden)
    with pytest.raises(bmt.exceptions.ArgumentError):
        project(result, cases['source'])
    assert result.to_dict() == before


def test_changed_reference_state_inventory_refuses_without_relabeling_original(
    interface, cases
):
    changed = msm.copy(cases['source'])
    changed.chemical_states.append_state()
    before = interface.to_dict()
    with pytest.raises(ValueError):
        project(interface, changed)
    assert interface.to_dict() == before


@pytest.mark.parametrize('missing', ['frame', 'state'])
def test_preserved_unknown_context_remains_readable_but_cannot_support_projection(
    interface, cases, missing
):
    record = interface.to_dict()
    metadata = record['analyses']['pisa']['topomt_observation']['metadata']
    metadata[
        'structure_indices' if missing == 'frame' else 'chemical_state_context'
    ] = None
    record['content_sha256'] = _digest(record)
    preserved = bmt.BindingInterface.from_dict(record, experimental=True)
    assert preserved.copy().to_dict() == record
    with pytest.raises(bmt.exceptions.ArgumentError):
        project(preserved, cases['source'])
    assert preserved.to_dict() == record


def test_source_mutation_during_owner_recovery_refuses_result(
    interface, cases, tmt, monkeypatch
):
    changed = msm.copy(cases['source'])
    load = tmt.load_provider_output

    def mutate_after_recovery(*args, **kwargs):
        output = load(*args, **kwargs)
        coordinates = msm.get(
            changed, selection=[4], structure_indices=[0], coordinates=True
        )
        values = puw.get_value(coordinates, to_unit='nm').copy()
        values[0, 0, 0] += 0.01
        msm.set(
            changed,
            selection=[4],
            structure_indices=[0],
            coordinates=puw.quantity(values, 'nm'),
        )
        return output

    monkeypatch.setattr(tmt, 'load_provider_output', mutate_after_recovery)
    before = interface.to_dict()
    with pytest.raises(bmt.exceptions.ArgumentError):
        project(interface, changed)
    assert interface.to_dict() == before


@pytest.mark.parametrize('kind', ['incomplete', 'duplicate', 'unsafe', 'corrupt'])
def test_actual_owner_recovery_rejects_derived_bad_archives_without_changing_original(
    interface, cases, kind, tmp_path
):
    original_bytes = (ORIGINAL / 'pisa.zip').read_bytes()
    path = tmp_path / f'{kind}.zip'
    if kind == 'corrupt':
        path.write_bytes(b'not a ZIP')
    else:
        with (
            zipfile.ZipFile(ORIGINAL / 'pisa.zip') as old,
            zipfile.ZipFile(path, 'w') as new,
        ):
            names = old.namelist()
            omitted = next(n for n in names if n != 'manifest.json')
            for name in names:
                if kind != 'incomplete' or name != omitted:
                    new.writestr(name, old.read(name))
            if kind == 'duplicate':
                with pytest.warns(UserWarning, match='Duplicate name'):
                    new.writestr('manifest.json', old.read('manifest.json'))
            elif kind == 'unsafe':
                new.writestr('../outside', b'unsafe control')
    before = interface.to_dict()
    with pytest.raises((ValueError, zipfile.BadZipFile)):
        interface.recover_topomt(cases['source'], bundle_file=path)
    assert interface.to_dict() == before
    assert (ORIGINAL / 'pisa.zip').read_bytes() == original_bytes


@pytest.mark.parametrize(
    'error',
    [
        MemoryError('controlled'),
        OSError(12, 'controlled'),
        KeyboardInterrupt('controlled'),
        SystemExit(4),
    ],
)
def test_original_owner_resource_and_process_errors_keep_identity(
    interface, cases, tmt, error, monkeypatch
):
    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(tmt, 'load_provider_output', fail)
    before = interface.to_dict()
    with pytest.raises(type(error)) as caught:
        project(interface, cases['source'])
    assert caught.value is error and interface.to_dict() == before


def test_new_receiver_versions_do_not_relabel_original_or_historical_records(
    interface, cases, tmp_path, monkeypatch
):
    footprint = project(interface, cases['source'])
    before = footprint.to_dict()
    old = interface.to_dict()
    interface.save(tmp_path / 'I.json')
    monkeypatch.setattr(bmt, '__version__', 'different-receiver')
    assert (
        bmt.BindingInterface.read(tmp_path / 'I.json', experimental=True).to_dict()
        == old
    )
    newer = intake(cases['positive'])
    assert (
        newer.analyses['pisa']['extensions']['bindingmt.original_intake'][
            'bindingmt_version'
        ]
        == 'different-receiver'
    )
    assert newer.topomt_analyses == interface.topomt_analyses
    check = footprint.check_context(cases['source'], bundle_file=ORIGINAL / 'pisa.zip')
    assert check['receiver']['bindingmt_version'] == 'different-receiver'
    assert footprint.to_dict() == before


def test_public_argument_admission_rejects_owner_context_override(interface, cases):
    with pytest.raises(UnknownArgumentError):
        project(interface, cases['source'], structure_index=7)


@pytest.mark.parametrize(
    'mutation',
    [
        'outcome',
        'run',
        'endpoint',
        'count',
        'measurement',
        'receiver',
        'asset',
        'version',
    ],
)
def test_resealed_contradictory_original_records_refuse(interface, mutation):
    record = interface.to_dict()
    analysis = record['analyses']['pisa']
    view = analysis['topomt_observation']
    if mutation == 'outcome':
        view['outcome'] = 'failed'
    elif mutation == 'run':
        analysis['reference']['run_id'] = 'a' * 64
    elif mutation == 'endpoint':
        view['findings'][0]['fields']['contacts']['other-bonds'][0][
            'source_atom_indices'
        ] = [400, 21]
    elif mutation == 'count':
        view['findings'][0]['fields']['contact_counts']['other-bonds'] = 0
    elif mutation == 'measurement':
        view['findings'][0]['measurements']['interface_area']['original_value'] += 1
    elif mutation == 'receiver':
        analysis['extensions']['bindingmt.original_intake']['original_execution'] = True
    elif mutation == 'asset':
        analysis['asset_refs'] = [{'uri': 'not-an-exact-reference'}]
    else:
        record['format_version'] = '0.5'
    record['content_sha256'] = _digest(record)
    original = deepcopy(record)
    with pytest.raises(bmt.exceptions.ArgumentError):
        bmt.BindingInterface.from_dict(record, experimental=True)
    assert record == original


def test_owned_endpoint_aggregation_refuses_unknown_or_many_to_one_labels():
    rows = [('g1', [{'source_atom_index': 2, 'contact_index': 0, 'endpoint_index': 0}])]
    for mapping in [{}, {2: 1, 7: 1}]:
        with pytest.raises(bmt.exceptions.ArgumentError):
            aggregate_endpoints(rows, mapping)
