"""Binding-specific preservation of a public PISA observation; no owner decoding."""

import hashlib
import json
import math
import re
from collections.abc import Mapping

import numpy as np
import pyunitwizard as puw
from pyunitwizard.record import QuantityRecord

from bindingmt._private.interface_records import _require
from bindingmt._private.serialization import json_copy

PROFILE = 'bindingmt.pisa_preserved_analysis@1'
GROUPS = ('h-bonds', 'salt-bridges', 'ss-bonds', 'cov-bonds', 'other-bonds')
INTERFACE_UNITS = {
    'interface_area': 'angstroms**2',
    'solvation_energy': 'kilocalories/mol',
    'stabilisation_energy': 'kilocalories/mol',
    'pvalue': 'dimensionless',
}
ASSEMBLY_UNITS = {
    'accessible_surface_area': 'angstroms**2',
    'buried_surface_area': 'angstroms**2',
    'dissociation_area': 'angstroms**2',
    'dissociation_energy': 'kilocalories/mol',
    'entropy_term': 'kilocalories/mol',
    'solvation_energy': 'kilocalories/mol',
}
MAX_CAPTURE_BYTES = 16 * 1024**2
MAX_RECORD_BYTES = 2 * 1024**2
MAX_FINDINGS = 256
MAX_CONTACTS = 65536


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _sha(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def _fields(value, names, field):
    _require(
        isinstance(value, Mapping) and set(value) == set(names),
        field,
        'Preserve the complete admitted PISA fields.',
    )


def measurement_field(scope, name):
    return f'bindingmt.topomt.pisa.{scope}.{name}'


def validate_receiver(value, operation):
    _fields(
        value,
        {'operation', 'bindingmt_version', 'topomt_version', 'original_execution'},
        'receiver',
    )
    _require(
        value['operation'] == operation
        and _text(value['bindingmt_version'])
        and _text(value['topomt_version'])
        and value['original_execution'] is False,
        'receiver',
        'Preserve the actual receiving operation/versions without originating-engine execution credit.',
    )


def measurements_to_wire(measurements, units, scope, run_id):
    _require(
        set(measurements) == set(units),
        'measurements',
        'Unsupported original PISA measurement inventory.',
    )
    result = {}
    for name, value in measurements.items():
        _require(
            value.run_id == run_id,
            'measurements',
            'Preserve original measurement run identity.',
        )
        result[name] = {
            'value': QuantityRecord.from_quantity(
                value.value, field=measurement_field(scope, name)
            ).to_dict(),
            'original_value': value.original_value,
            'original_unit': value.original_unit,
            'source_field': value.source_field,
            'source_artifact': value.source_artifact,
            'run_id': value.run_id,
            'definition': value.definition,
            'issue_url': value.issue_url,
            'status': value.status,
        }
    return result


def _validate_measurements(values, units, scope, run_id):
    _fields(values, units, 'measurements')
    for name, value in values.items():
        _fields(
            value,
            {
                'value',
                'original_value',
                'original_unit',
                'source_field',
                'source_artifact',
                'run_id',
                'definition',
                'issue_url',
                'status',
            },
            'measurements',
        )
        unit = units[name]
        _require(
            value['original_unit'] == unit
            and value['run_id'] == run_id
            and value['status'] == 'external_only',
            'measurements',
            'Preserve original PISA units, identity and external-only authority.',
        )
        _require(
            all(
                _text(value[k])
                for k in ('source_field', 'source_artifact', 'definition')
            )
            and (value['issue_url'] is None or _text(value['issue_url'])),
            'measurements',
            'Preserve original field, artifact and definition references.',
        )
        original = value['original_value']
        _require(
            type(original) in (int, float) and math.isfinite(original),
            'measurements',
            'Use the admitted finite scalar PISA measures.',
        )
        _require(
            (unit != 'angstroms**2' or original >= 0)
            and (name != 'pvalue' or 0 <= original <= 1),
            'measurements',
            'Preserve nonnegative original areas and the reported probability range.',
        )
        quantity = QuantityRecord.from_dict(value['value']).to_quantity(
            field=measurement_field(scope, name),
            unit=unit,
            dimensionality=puw.get_dimensionality(puw.unit(unit, form='pint')),
            form='pint',
        )
        number = np.asarray(puw.get_value(quantity, to_unit=unit))
        _require(
            number.shape == ()
            and math.isclose(
                float(number), original, rel_tol=32 * np.finfo(float).eps, abs_tol=0
            ),
            'measurements',
            'The quantity and original reported value disagree.',
        )


def snapshot(provider_result):
    """Called only after the optional owner's public types have been admitted."""
    run = provider_result.run if hasattr(provider_result, 'run') else provider_result
    _require(
        sum(len(data) for _, data in run.artifacts) <= MAX_CAPTURE_BYTES,
        'provider_result',
        'Original capture exceeds the admitted 16 MiB input bound.',
    )
    decoded = hasattr(provider_result, 'interfaces')
    if decoded:
        _require(
            len(provider_result.records)
            == len(provider_result.interfaces)
            <= MAX_FINDINGS,
            'provider_result',
            'Use at most 256 PISA interface records.',
        )
        count = sum(
            len(rows)
            for record in provider_result.interfaces
            for rows in record.fields.get('contacts', {}).values()
        )
        _require(
            count <= MAX_CONTACTS,
            'provider_result',
            'Original contacts exceed the admitted bound.',
        )
        findings = []
        for record in provider_result.interfaces:
            _require(
                record.run == run and not record.geometries,
                'provider_result',
                'Use the admitted original PISA run without an unqualified geometry profile.',
            )
            fields = record.fields
            fields.pop('external_measurements', None)
            findings.append(
                {
                    'source_id': record.source_id,
                    'atom_role': record.atom_role,
                    'atom_indices': None
                    if record.atom_indices is None
                    else list(record.atom_indices),
                    'atom_labels': None
                    if record.atom_labels is None
                    else list(record.atom_labels),
                    'parent_source_ids': list(record.parent_source_ids),
                    'geometry_availability': 'not_supplied',
                    'fields': fields,
                    'measurements': measurements_to_wire(
                        record.measurements,
                        INTERFACE_UNITS,
                        f'interface.{record.source_id}',
                        run.run_id,
                    ),
                }
            )
        outcome = 'completed'
        assembly = measurements_to_wire(
            provider_result.assembly_measurements,
            ASSEMBLY_UNITS,
            'assembly',
            run.run_id,
        )
    else:
        findings, assembly = None, None
        outcome = (
            'failed'
            if run.metadata.get('execution_status') == 'failed'
            and run.metadata.get('failure_stage') == 'analysis'
            else 'unknown'
        )
    result = json_copy(
        {
            'profile': PROFILE,
            'provider': run.provider,
            'backend': run.backend,
            'run_id': run.run_id,
            'metadata': run.metadata,
            'artifact_checksums': {
                name: hashlib.sha256(data).hexdigest() for name, data in run.artifacts
            },
            'attribution': run.attribution,
            'outcome': outcome,
            'assembly_policy': provider_result.assembly_policy
            if decoded
            else run.metadata.get('assembly_policy'),
            'assembly_measurements': assembly,
            'findings': findings,
        },
        'provider_result',
    )
    validate_observation(result)
    return result


def validate_observation(value):
    """Validate only this owned, detached wire profile; import no optional owner."""
    value = json_copy(value, 'topomt_observation')
    _fields(
        value,
        {
            'profile',
            'provider',
            'backend',
            'run_id',
            'metadata',
            'artifact_checksums',
            'attribution',
            'outcome',
            'assembly_policy',
            'assembly_measurements',
            'findings',
        },
        'topomt_observation',
    )
    _require(
        len(json.dumps(value, allow_nan=False).encode()) <= MAX_RECORD_BYTES,
        'topomt_observation',
        'Use at most 2 MiB of preserved observation.',
    )
    _require(
        value['profile'] == PROFILE
        and value['provider'] == 'pisa'
        and value['backend'] == 'cli'
        and _sha(value['run_id']),
        'topomt_observation',
        'Use the admitted PISA run/profile identity.',
    )
    _require(
        isinstance(value['metadata'], Mapping)
        and value['metadata'].get('adapter_protocol') == 'topomt.original-local@1',
        'topomt_observation',
        'Preserve the original local PISA protocol context.',
    )
    selected = value['metadata'].get('selected_atom_indices')
    _require(
        isinstance(selected, list)
        and selected
        and all(type(i) is int and i >= 0 for i in selected)
        and len(selected) == len(set(selected)),
        'topomt_observation',
        'Preserve the original selected input axis.',
    )
    checksums = value['artifact_checksums']
    _require(
        isinstance(checksums, Mapping)
        and bool(checksums)
        and all(_text(k) and _sha(v) for k, v in checksums.items()),
        'topomt_observation',
        'Preserve original artifact checksum references.',
    )
    _require(
        value['attribution'] is None or isinstance(value['attribution'], Mapping),
        'attribution',
        'Preserve the original attribution or its explicit absence.',
    )
    outcome = value['outcome']
    _require(
        outcome in ('completed', 'failed', 'unknown')
        and value['assembly_policy'] == 'supplied_as_is',
        'topomt_observation',
        'Preserve the admitted assembly and original outcome.',
    )
    if outcome != 'completed':
        _require(
            value['findings'] is None and value['assembly_measurements'] is None,
            'topomt_observation',
            'Unrecovered/failed evidence has unavailable findings, not completed zero.',
        )
        if outcome == 'failed':
            _require(
                value['metadata'].get('execution_status') == 'failed'
                and value['metadata'].get('failure_stage') == 'analysis',
                'topomt_observation',
                'Original analysis failure needs its confirmed owner metadata.',
            )
        return
    _require(
        value['metadata'].get('execution_status') == 'completed',
        'topomt_observation',
        'Completed recovery requires the original completed execution context.',
    )
    _validate_measurements(
        value['assembly_measurements'], ASSEMBLY_UNITS, 'assembly', value['run_id']
    )
    findings = value['findings']
    _require(
        isinstance(findings, list) and len(findings) <= MAX_FINDINGS,
        'topomt_observation',
        'Preserve the bounded ordered finding list, including completed [].',
    )
    identities, total = [], 0
    for finding in findings:
        _fields(
            finding,
            {
                'source_id',
                'atom_role',
                'atom_indices',
                'atom_labels',
                'parent_source_ids',
                'geometry_availability',
                'fields',
                'measurements',
            },
            'findings',
        )
        identity = finding['source_id']
        _require(
            _text(identity) and identity not in identities,
            'findings',
            'Preserve distinct original finding identities.',
        )
        identities.append(identity)
        _require(
            finding['geometry_availability'] == 'not_supplied'
            and finding['parent_source_ids'] == []
            and (
                finding['atom_labels'] is None
                or isinstance(finding['atom_labels'], list)
            ),
            'findings',
            'Preserve the admitted PISA geometry/parent availability.',
        )
        indices = finding['atom_indices']
        _require(
            indices is None
            or isinstance(indices, list)
            and all(type(i) is int and i >= 0 for i in indices)
            and indices == sorted(set(indices)),
            'findings',
            'Preserve exact canonical source atom indices or unavailable mapping.',
        )
        fields = finding['fields']
        _require(
            isinstance(fields, Mapping)
            and fields.get('source_id') == identity
            and isinstance(fields.get('participants'), list)
            and len(fields['participants']) == 2,
            'findings',
            'Preserve the finding and its two original participants.',
        )
        chains = [
            p.get('chain_id') if isinstance(p, Mapping) else None
            for p in fields['participants']
        ]
        _require(
            all(_text(c) for c in chains) and len(set(chains)) == 2,
            'findings',
            'Preserve distinct reported participant chains.',
        )
        _fields(fields.get('contacts'), GROUPS, 'contacts')
        _fields(fields.get('contact_counts'), GROUPS, 'contacts')
        members = set()
        for group, rows in fields['contacts'].items():
            _require(
                isinstance(rows, list)
                and type(fields['contact_counts'][group]) is int
                and fields['contact_counts'][group] == len(rows),
                'contacts',
                'Preserve every original contact row and its count.',
            )
            total += len(rows)
            _require(
                total <= MAX_CONTACTS,
                'contacts',
                'Original contacts exceed the admitted bound.',
            )
            for row in rows:
                _require(
                    isinstance(row, Mapping)
                    and isinstance(row.get('original'), Mapping)
                    and isinstance(row.get('distance'), Mapping),
                    'contacts',
                    'Preserve original rows and unmodified distance packets.',
                )
                pair = row.get('source_atom_indices')
                if finding['atom_indices'] is None:
                    _require(
                        pair is None,
                        'contacts',
                        'An unavailable map cannot silently contain partially mapped endpoints.',
                    )
                else:
                    _require(
                        isinstance(pair, list)
                        and len(pair) == 2
                        and all(type(i) is int and i >= 0 for i in pair),
                        'contacts',
                        'Preserve exact mapped source endpoints.',
                    )
                    members.update(pair)
        expected_role = (
            'unavailable_contact_atom_map'
            if finding['atom_indices'] is None
            else 'reported_contact_atoms'
        )
        _require(
            finding['atom_role'] == expected_role
            and (
                finding['atom_indices'] is None
                or finding['atom_indices'] == sorted(members)
            ),
            'findings',
            'Original contact membership, role and complete row witnesses disagree.',
        )
        _require(
            members <= set(selected),
            'findings',
            'Reported endpoints must belong to the original submitted source selection.',
        )
        _validate_measurements(
            finding['measurements'],
            INTERFACE_UNITS,
            f'interface.{identity}',
            value['run_id'],
        )


def finding_refs(value):
    return (
        [
            {'run_id': value['run_id'], 'source_id': finding['source_id']}
            for finding in value['findings']
        ]
        if value['outcome'] == 'completed'
        else None
    )
