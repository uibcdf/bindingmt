"""Independent immutable footprints of native or reported contact observations."""

import json
from dataclasses import dataclass

from argdigest import arg_digest

from bindingmt._private.interface_contacts import (
    footprint_payload,
    validate_calculation,
)
from bindingmt._private.interface_records import _digest, _InterfaceRecordIO, _require
from bindingmt._private.serialization import json_copy
from bindingmt._private.smonitor import signal


@dataclass(frozen=True, init=False, repr=False)
class ContactFootprint(_InterfaceRecordIO):
    """Detached source-ordered members and all atomic witnesses, not a binding site."""

    _json: str

    @arg_digest(config='bindingmt._argdigest')
    def __init__(self, record, *, experimental=False):
        _require(
            experimental is True,
            'experimental',
            'Opt in to experimental footprint 0.1/0.2.',
        )
        version = record.get('format_version')
        body_field = 'reported_projection' if version == '0.2' else 'calculation'
        required = {
            'format',
            'format_version',
            'experimental',
            'record_id',
            'content_sha256',
            'parent_record_id',
            body_field,
            'footprint',
        }
        _require(
            set(record) == required
            and record['format'] == 'bindingmt.contact_footprint'
            and version in ('0.1', '0.2')
            and record['experimental'] is True,
            'record',
            'Use a complete experimental contact footprint 0.1/0.2 record.',
        )
        _require(
            all(
                isinstance(record[k], str) and record[k]
                for k in ('record_id', 'parent_record_id')
            ),
            'record_id',
            'Preserve footprint and parent identities.',
        )
        _require(
            record['record_id'] != record['parent_record_id'],
            'parent_record_id',
            'A footprint cannot name itself as its parent.',
        )
        snapshot = json_copy(record, 'record')
        if version == '0.1':
            validate_calculation(snapshot['calculation'])
            expected = (
                footprint_payload(
                    snapshot['calculation'], snapshot['footprint'].get('participant')
                )
                if isinstance(snapshot['footprint'], dict)
                else None
            )
        else:
            from bindingmt._private.interface_reported_contacts import (
                validate_projection,
            )

            expected = validate_projection(snapshot['reported_projection'])
            _require(
                snapshot['reported_projection']['parent_ref']['record_id']
                == snapshot['parent_record_id'],
                'parent_record_id',
                'Preserve the reported projection parent identity.',
            )
        payload = snapshot['footprint']
        _require(
            isinstance(payload, dict)
            and json.dumps(payload, sort_keys=True)
            == json.dumps(expected, sort_keys=True),
            'footprint',
            'Footprint members, order/maps and contributors disagree with the calculation.',
        )
        _require(
            snapshot['content_sha256'] == _digest(snapshot),
            'content_sha256',
            'Footprint content digest mismatch.',
        )
        object.__setattr__(self, '_json', json.dumps(snapshot, allow_nan=False))

    @property
    def source_atom_indices(self):
        return self.to_dict()['footprint']['source_atom_indices']

    @property
    def local_atom_indices(self):
        return self.to_dict()['footprint']['local_atom_indices']

    def summary(self):
        record = self.to_dict()
        if record['format_version'] == '0.2':
            projection = record['reported_projection']
            return {
                **record['footprint'],
                'record_id': record['record_id'],
                'parent_record_id': record['parent_record_id'],
                'analysis_ref': projection['analysis_ref'],
                'finding_id': projection['finding_id'],
                'definition': projection['definition'],
                'context': projection['source_context'],
                'current_assessment': 'not_performed',
                'asset_access_on_read': 'not_performed',
            }
        return {
            **record['footprint'],
            'record_id': record['record_id'],
            'parent_record_id': record['parent_record_id'],
            'analysis_id': record['calculation']['analysis_id'],
            'context': record['calculation']['context'],
            'current_assessment': 'not_performed',
        }

    @signal(tags=['api', 'interface', 'context'])
    @arg_digest(config='bindingmt._argdigest')
    def check_context(self, molecular_system, *, bundle_file):
        """Explicitly reassess a reported footprint without mutating its history."""
        from bindingmt.provider_interfaces import _check_reported_context

        return _check_reported_context(self, molecular_system, bundle_file=bundle_file)
