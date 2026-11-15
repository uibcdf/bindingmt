"""Exercise the copyable example itself, including isolated historical reading."""

import hashlib
import json
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import pyunitwizard as puw

import bindingmt as bmt


def example_path():
    # Installed tests receive these exact files through the shared SDK namespace.
    import bmt_test_support

    staged = Path(bmt_test_support.__path__[0]) / 'site_interface/walkthrough.py'
    if staged.is_file():
        return staged
    return (
        Path(__file__).resolve().parents[1] / 'examples/site_interface/walkthrough.py'
    )


@pytest.fixture
def example():
    return runpy.run_path(str(example_path()))


@pytest.mark.parametrize('frame', (0, 1, 2))
@pytest.mark.parametrize(
    'form,units',
    [('pint', ['nm', 'ps']), ('pint', ['pm', 'fs']), ('unyt', ['pm', 'fs'])],
)
def test_actual_walkthrough_observations_and_policy(
    example, tmp_path, frame, form, units
):
    data = example_path().parent
    source_before = (data / 'source.h5msm').read_bytes()
    expected = json.loads((data / 'case.json').read_text())['expected'][frame]
    with puw.context(default_form=form, default_parser='pint', standard_units=units):
        policy = puw.configure.report()
        output = tmp_path / 'new-run'
        result = example['run_example'](output, structure_index=frame)
        assert puw.configure.report() == policy
        region = result['region']['projection']
        assert region['exported_target_atom_indices'] == expected['region']
        assert region['axes']['source_structure_index'] == frame
        assert region['axes']['chemical_state_index'] == expected['state']
        assert (
            result['contacts'].participants['partner']['source_atom_indices']
            == expected['partner']
        )
        assert (
            result['contacts'].contact_calculation['local_pairs'] == expected['pairs']
        )
        assert result['contacts'].analyses['contacts']['outcome'] == 'completed'
        assert result['footprint'].source_atom_indices == expected['footprint']
        for key in ('contacts', 'burial'):
            context = result[key].to_dict()['contexts']['calculation']
            assert context['source_structure_index'] == frame
            assert context['structure_chemical_state_index'] == expected['state']
        for area in result['burial'].burial_areas['sides'].values():
            value = float(puw.get_value(area, to_unit='nm**2'))
            assert value > 0 if expected['burial'] == 'positive' else value == 0
        assert result['burial'].summary()['native_surface']['burial_usable'] is True
        restored = example['read_results'](output)
        for key in ('site', 'contacts', 'burial', 'footprint'):
            assert restored[key].to_dict() == result[key].to_dict()
        assert restored['region'] == result['region']
        assert puw.configure.report() == policy
    assert (data / 'source.h5msm').read_bytes() == source_before
    assert (output / 'source.h5msm').read_bytes() == source_before


@pytest.mark.parametrize('frame', (1, 2))
def test_copied_cli_and_cold_reader_need_no_original_source(tmp_path, frame):
    copied = tmp_path / 'copied'
    copied.mkdir()
    for name in ('walkthrough.py', 'case.json', 'source.h5msm'):
        shutil.copyfile(example_path().parent / name, copied / name)
    output = tmp_path / 'results'
    run = subprocess.run(
        [
            sys.executable,
            '-I',
            '-B',
            str(copied / 'walkthrough.py'),
            str(output),
            '--structure-index',
            str(frame),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )
    assert 'Candidate region: [2]' in run.stdout
    footprint = '[2]' if frame == 1 else '[]'
    assert f'Observed region footprint: {footprint}' in run.stdout
    assert 'nm**2' in run.stdout
    before = {p.name: p.read_bytes() for p in output.iterdir()}
    repeated = subprocess.run(
        [sys.executable, '-I', '-B', str(copied / 'walkthrough.py'), str(output)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert repeated.returncode == 2
    assert 'fresh output directory' in repeated.stderr
    assert {p.name: p.read_bytes() for p in output.iterdir()} == before
    for path in (
        copied / 'case.json',
        copied / 'source.h5msm',
        output / 'source.h5msm',
    ):
        path.unlink()
    reader = tmp_path / 'cold.py'
    reader.write_text("""import runpy
import sys
import molsysmt as msm
import bindingmt as bmt
def forbidden(*args, **kwargs):
    raise AssertionError('Historical read attempted molecular work')
for name in ('get', 'set', 'select', 'convert', 'extract', 'copy'):
    setattr(msm, name, forbidden)
msm.h5msm.read = forbidden
msm.structure.get_contacts = forbidden
msm.physchem.get_sasa = forbidden
bmt.calculate_interface_contacts = forbidden
bmt.calculate_interface_burial = forbidden
sys.argv = [sys.argv[1], sys.argv[2], '--read']
runpy.run_path(sys.argv[0], run_name='__main__')
assert 'vina' not in sys.modules
assert 'molsysviewer' not in sys.modules
""")
    cold = subprocess.run(
        [
            sys.executable,
            '-I',
            '-B',
            str(reader),
            str(copied / 'walkthrough.py'),
            str(output),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )
    assert 'Saved history; molecular source not accessed.' in cold.stdout
    assert f'Observed region footprint: {footprint}' in cold.stdout
    assert {p.name: p.read_bytes() for p in output.iterdir()} == {
        name: value for name, value in before.items() if name != 'source.h5msm'
    }


def test_saved_site_keeps_its_captured_structure(example, tmp_path):
    import molsysmt as msm

    result = example['run_example'](tmp_path / 'frame-1', structure_index=1)
    source = msm.h5msm.read(str(example_path().parent / 'source.h5msm'))
    original = result['site'].to_dict()
    with pytest.raises(bmt.exceptions.ArgumentError):
        result['site'].resolve_region(source, structure_index=2)
    assert result['site'].to_dict() == original


@pytest.mark.parametrize('mode', ('missing_receipt', 'changed_file'))
def test_reader_refuses_incomplete_or_changed_results(example, tmp_path, mode):
    output = tmp_path / 'results'
    example['run_example'](output)
    if mode == 'missing_receipt':
        (output / 'run.json').unlink()
        with pytest.raises(FileNotFoundError):
            example['read_results'](output)
    else:
        with (output / 'contacts.json').open('a') as stream:
            stream.write('\n')
        with pytest.raises(ValueError, match='differs from its run receipt'):
            example['read_results'](output)


def test_required_calculation_failure_is_preserved_without_completed_receipt(
    example, tmp_path, monkeypatch
):
    error = OSError(12, 'controlled required child resource failure')

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(bmt, 'calculate_interface_burial', fail)
    output = tmp_path / 'failed-run'
    with pytest.raises(OSError) as caught:
        example['run_example'](output)
    assert caught.value is error
    assert not (output / 'run.json').exists()
    assert not list(output.glob('*.json'))


def test_save_failure_keeps_partial_records_without_completed_receipt(
    example, tmp_path, monkeypatch
):
    error = OSError(5, 'controlled footprint save failure')

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(bmt.ContactFootprint, 'save', fail)
    output = tmp_path / 'partial-run'
    with pytest.raises(OSError) as caught:
        example['run_example'](output)
    assert caught.value is error
    assert (output / 'site.json').is_file()
    assert (output / 'contacts.json').is_file()
    assert (output / 'burial.json').is_file()
    assert not (output / 'run.json').exists()


def test_frozen_input_and_authoring_identity():
    import bmt_test_support.build_site_interface_input as author

    data = example_path().parent
    case = json.loads((data / 'case.json').read_text())
    payload = (data / 'source.h5msm').read_bytes()
    assert hashlib.sha256(payload).hexdigest() == case['source']['sha256']
    assert len(payload) == case['source']['bytes']
    assert (
        hashlib.sha256(Path(author.__file__).read_bytes()).hexdigest()
        == case['producer']['authoring_program_sha256']
    )
