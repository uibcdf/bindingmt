"""Exact offered PISA inputs for source and officially staged test support."""

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / 'tests/data/interface_topomt_pisa'
ORIGINAL = DATA / '3a8ee979'
NATIVE = DATA / 'native_6deb0200'
RECORDS = DATA.parent / 'interface_topomt_records/20261115'


def verify():
    manifest = json.loads((DATA / 'receiving.json').read_text())
    for name, expected in manifest['files'].items():
        assert hashlib.sha256((ORIGINAL / name).read_bytes()).hexdigest() == expected
    native = json.loads((NATIVE / 'manifest.json').read_text())
    for value in native['sources'].values():
        payload = (NATIVE / value['file']).read_bytes()
        assert len(payload) == value['bytes']
        assert hashlib.sha256(payload).hexdigest() == value['sha256']
    return manifest


def source(*, empty=False):
    import molsysmt as msm

    verify()
    return msm.h5msm.read(str(NATIVE / ('empty.h5msm' if empty else 'positive.h5msm')))


def verify_records():
    manifest = json.loads((RECORDS / 'manifest.json').read_text())
    for name, expected in manifest['files'].items():
        assert hashlib.sha256((RECORDS / name).read_bytes()).hexdigest() == expected
    return manifest
