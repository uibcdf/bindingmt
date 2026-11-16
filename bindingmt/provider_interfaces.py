"""Explicit optional intake of original interface evidence through its owner."""

import hashlib
import json
import os
from collections.abc import Mapping
from pathlib import Path
from uuid import uuid4

from depdigest import dep_digest

from bindingmt._private.interface_records import _digest, _require
from bindingmt._private.interface_topomt import (
    MAX_CAPTURE_BYTES,
    finding_refs,
    snapshot,
)
from bindingmt._private.serialization import json_copy


def _receiver(operation):
    import topomt as tmt

    from bindingmt import __version__

    return {
        'operation': operation,
        'bindingmt_version': __version__,
        'topomt_version': tmt.__version__,
        'original_execution': False,
    }


@dep_digest('topomt')
def _from_topomt(
    provider_result,
    *,
    question,
    participants,
    grouping,
    analysis_name,
    artifact_ref,
    contexts,
    site_refs,
):
    """Consume public owner objects; neither decode files nor execute an engine."""
    from topomt.provider_output import ProviderRun
    from topomt.third_party.pisa.output import PisaOutput

    from bindingmt.core.interface import BindingInterface, characterize_interface

    _require(
        isinstance(provider_result, (PisaOutput, ProviderRun)),
        'provider_result',
        'Supply a public PisaOutput or ProviderRun; no foreign/duck-typed output is admitted.',
    )
    _require(
        isinstance(analysis_name, str) and bool(analysis_name.strip()),
        'analysis_name',
        'Name the preserved analysis.',
    )
    if artifact_ref is not None:
        from bindingmt._private.interface_topomt import _sha

        _require(
            set(artifact_ref) == {'uri', 'sha256'}
            and isinstance(artifact_ref['uri'], str)
            and bool(artifact_ref['uri'].strip())
            and _sha(artifact_ref['sha256']),
            'artifact_ref',
            'Declare an exact URI/SHA-256 reference; it is not accessed here.',
        )
    observation = snapshot(provider_result)
    base = characterize_interface(
        question=question,
        participants=participants,
        grouping=grouping,
        contexts=contexts,
        site_refs=site_refs,
    )
    record = base.to_dict()
    original = {
        key: record[key] for key in ('question', 'participants', 'grouping', 'contexts')
    }
    record['format_version'] = '0.6'
    record['analyses'] = {
        analysis_name: {
            'owner': 'pisa',
            'reference': {
                key: observation[key] for key in ('provider', 'backend', 'run_id')
            },
            'authority': 'topomt_preserved',
            'outcome': observation['outcome'],
            'finding_refs': finding_refs(observation),
            'original_claim': original,
            'asset_refs': [] if artifact_ref is None else [artifact_ref],
            'topomt_observation': observation,
            'extensions': {
                'bindingmt.original_intake': _receiver('BindingInterface.from_topomt')
            },
        }
    }
    if observation['outcome'] == 'failed':
        record['analyses'][analysis_name]['failure'] = {
            'stage': 'analysis',
            'reference': {'run_id': observation['run_id']},
            'artifact_refs': sorted(observation['artifact_checksums']),
        }
    record['content_sha256'] = _digest(record)
    return BindingInterface(record, experimental=True)


def _same_original_facts(left, right):
    """Compare exact owner facts after each typed quantity validates independently.

    QuantityRecord wire presentation may differ across unit backends/policies;
    each value is already checked against its same original scalar and unit.
    No opaque provenance, version or source field is filtered.
    """

    def facts(value):
        result = json_copy(value, 'topomt_observation')
        measures = [result['assembly_measurements']]
        measures.extend(f['measurements'] for f in result['findings'] or [])
        for collection in measures:
            if collection is not None:
                for measurement in collection.values():
                    del measurement['value']
        return result

    return json.dumps(facts(left), sort_keys=True) == json.dumps(
        facts(right), sort_keys=True
    )


@dep_digest('topomt')
def _recover(observation, molecular_system, bundle_file, artifact_refs=()):
    import molsysmt as msm
    import topomt as tmt
    from molsysmt.native import MolSys

    from bindingmt._private.interface_context import _geometry_digest
    from bindingmt._private.interface_reported_contacts import (
        CHECKS,
        MAX_SOURCE_ATOMS,
        PRECONDITIONS,
    )

    _require(
        observation['outcome'] == 'completed',
        'analysis_name',
        'Only a completed original analysis has this source-bound recovery profile.',
    )
    _require(
        isinstance(molecular_system, MolSys),
        'molecular_system',
        'Supply the complete original native MolSys.',
    )
    count, frames = msm.get(molecular_system, n_atoms=True, n_structures=True)
    frame = observation['metadata'].get('structure_indices')
    state_context = observation['metadata'].get('chemical_state_context')
    _require(
        isinstance(state_context, Mapping),
        'molecular_system',
        'Original reference-state correspondence is unavailable.',
    )
    state = state_context.get('index')
    _require(
        0 < count <= MAX_SOURCE_ATOMS
        and type(frame) is int
        and 0 <= frame < frames
        and type(state) is int
        and state >= 0
        and observation['metadata'].get('chemical_state') == 'reference',
        'molecular_system',
        'Use the original bounded source, explicit saved frame and reference chemical state.',
    )
    _require(
        isinstance(bundle_file, (str, os.PathLike)),
        'bundle_file',
        'Supply an explicit local original bundle path; no automatic asset lookup.',
    )
    path = Path(bundle_file)
    _require(
        path.stat().st_size <= MAX_CAPTURE_BYTES,
        'bundle_file',
        'The bundle exceeds this admitted input bound.',
    )
    payload = path.read_bytes()
    _require(
        len(payload) <= MAX_CAPTURE_BYTES,
        'bundle_file',
        'The bundle changed beyond this admitted input bound.',
    )
    digest = hashlib.sha256(payload).hexdigest()
    for reference in artifact_refs:
        _require(
            reference.get('sha256') == digest,
            'bundle_file',
            'The supplied bytes differ from the preserved exact asset reference.',
        )
    before = _geometry_digest(msm, molecular_system, list(range(count)), frame)
    output = tmt.load_provider_output(str(path), molecular_system=molecular_system)
    recovered = snapshot(output)
    _require(
        _same_original_facts(observation, recovered),
        'bundle_file',
        'The recovered original analysis differs from the preserved evidence.',
    )
    _require(
        path.read_bytes() == payload
        and _geometry_digest(msm, molecular_system, list(range(count)), frame)
        == before,
        'molecular_system',
        'The original source or artifact changed during recovery.',
    )
    return output, {
        'scope': 'molecular_correspondence',
        'context_policy': 'require_match',
        'source_atom_count': int(count),
        'source_structure_index': frame,
        'chemical_state_index': state,
        'geometry_sha256_nm_float64': before,
        'bundle_sha256': digest,
        'performed_checks': CHECKS,
        'caller_preconditions': PRECONDITIONS,
        'assessment_kind': 'historical_check_receipt',
    }


def _analysis(interface, analysis_name):
    analyses = interface.topomt_analyses
    _require(
        analysis_name in analyses,
        'analysis_name',
        'Select an actual preserved original PISA analysis.',
    )
    return analyses[analysis_name]


def _recover_topomt(interface, molecular_system, *, analysis_name, bundle_file):
    observation = _analysis(interface, analysis_name)
    refs = interface.analyses[analysis_name].get('asset_refs', [])
    return _recover(observation, molecular_system, bundle_file, refs)[0]


def _bindings(
    molecular_system, observation, finding_id, participants, grouping, context
):
    import molsysmt as msm

    from bindingmt._private.interface_reported_contacts import finding

    original = finding(observation, finding_id)
    _require(
        len(participants) == 2
        and len(grouping) == 2
        and all(len(g) == 1 for g in grouping),
        'grouping',
        'This first PISA molecular profile binds two singleton participants.',
    )
    chains = [p['chain_id'] for p in original['fields']['participants']]
    atom_chains = msm.get(
        molecular_system, element='atom', selection='all', chain_name=True
    )
    _require(
        isinstance(atom_chains, (list, tuple))
        and len(atom_chains) == context['source_atom_count'],
        'molecular_system',
        'Original atom/chain correspondence is unavailable.',
    )
    axes = {
        chain: [i for i, value in enumerate(atom_chains) if value == chain]
        for chain in chains
    }
    bindings, ids = {}, []
    for group in grouping:
        name = group[0]
        reference = participants[name]
        indices = reference.get('source_atom_indices')
        source_id = reference.get('source_id')
        _require(
            isinstance(source_id, str)
            and source_id.strip()
            and isinstance(indices, list)
            and indices
            and all(type(i) is int for i in indices)
            and len(indices) == len(set(indices)),
            'participants',
            'Preserve explicit complete original chain axes and a declared source identity.',
        )
        matches = [
            chain for chain, axis in axes.items() if axis and set(axis) == set(indices)
        ]
        _require(
            len(matches) == 1,
            'participants',
            'Each declared participant must match one complete original PISA chain.',
        )
        bindings[name] = {
            'owner_chain_id': matches[0],
            'source_atom_indices': sorted(indices),
        }
        ids.append(source_id)
    _require(
        ids[0] == ids[1] and len({b['owner_chain_id'] for b in bindings.values()}) == 2,
        'participants',
        'Preserve a shared original source and distinct original chains.',
    )
    return bindings, ids[0]


def _project_reported(
    interface,
    molecular_system,
    *,
    bundle_file,
    analysis_name,
    finding_id,
    participant,
    contact_groups,
    assurance_scope,
    context_policy,
):
    from bindingmt._private.interface_reported_contacts import (
        DEFINITION,
        distance_records,
        finding,
        selected_groups,
        validate_projection,
    )
    from bindingmt._private.smonitor import CapabilityMismatchError
    from bindingmt.core.contact_footprint import ContactFootprint

    if assurance_scope != 'molecular_correspondence':
        raise CapabilityMismatchError(
            capability=assurance_scope, engine='bindingmt.reported_contacts'
        )
    _require(
        context_policy == 'require_match',
        'context_policy',
        'This molecular projection requires a fresh original-source match.',
    )
    groups = selected_groups(contact_groups)
    observation = _analysis(interface, analysis_name)
    original = finding(observation, finding_id)
    _require(
        original['atom_indices'] is not None,
        'mapping',
        'Original endpoint correspondence is unavailable, not evaluated empty.',
    )
    refs = interface.analyses[analysis_name].get('asset_refs', [])
    _, context = _recover(observation, molecular_system, bundle_file, refs)
    bindings, source_id = _bindings(
        molecular_system,
        observation,
        finding_id,
        interface.participants,
        interface.grouping,
        context,
    )
    _require(
        participant in bindings,
        'participant',
        'Select one of the admitted declared participants.',
    )
    projection = {
        'receiver': _receiver('BindingInterface.project_reported_footprint'),
        'definition': DEFINITION,
        'analysis_name': analysis_name,
        'analysis_ref': {k: observation[k] for k in ('provider', 'backend', 'run_id')},
        'parent_ref': {
            'record_id': interface.record_id,
            'content_sha256': interface.to_dict()['content_sha256'],
        },
        'observation': observation,
        'finding_id': finding_id,
        'source_id': source_id,
        'participant': participant,
        'participant_bindings': bindings,
        'contact_groups': groups,
        'local_atom_map': [[i, i] for i in range(context['source_atom_count'])],
        'source_context': context,
        'reported_distances': distance_records(original, groups),
    }
    _, final = _recover(observation, molecular_system, bundle_file, refs)
    _require(
        final == context,
        'molecular_system',
        'Original source/artifact context changed during projection.',
    )
    footprint = validate_projection(projection)
    record = {
        'format': 'bindingmt.contact_footprint',
        'format_version': '0.2',
        'experimental': True,
        'record_id': str(uuid4()),
        'parent_record_id': interface.record_id,
        'reported_projection': projection,
        'footprint': footprint,
    }
    record['content_sha256'] = _digest(record)
    return ContactFootprint(record, experimental=True)


def _check_reported_context(footprint, molecular_system, *, bundle_file):
    from bindingmt._private.interface_reported_contacts import validate_projection

    record = footprint.to_dict()
    _require(
        record['format_version'] == '0.2',
        'footprint',
        'This explicit context checker admits the reported PISA footprint profile only.',
    )
    projection = record['reported_projection']
    validate_projection(projection)
    observation = projection['observation']
    _, context = _recover(observation, molecular_system, bundle_file)
    grouping = [[name] for name in projection['participant_bindings']]
    participants = {
        name: {
            'source_id': projection['source_id'],
            'source_atom_indices': b['source_atom_indices'],
        }
        for name, b in projection['participant_bindings'].items()
    }
    bindings, _ = _bindings(
        molecular_system,
        observation,
        projection['finding_id'],
        participants,
        grouping,
        context,
    )
    _require(
        bindings == projection['participant_bindings']
        and context == projection['source_context'],
        'molecular_system',
        'The proposed current source/artifact does not match this historical footprint context.',
    )
    _, final = _recover(observation, molecular_system, bundle_file)
    _require(
        final == context,
        'molecular_system',
        'The proposed current source changed during checking.',
    )
    return {
        'footprint_ref': {
            'record_id': record['record_id'],
            'content_sha256': record['content_sha256'],
        },
        'source_context': context,
        'result': 'match',
        'receiver': _receiver('ContactFootprint.check_context'),
    }
