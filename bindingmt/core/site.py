"""Immutable, coordinate-free binding-site hypotheses and source annotations."""

import json
from dataclasses import dataclass
from uuid import uuid4

from argdigest import arg_digest

from bindingmt._private.interface_records import _digest, _InterfaceRecordIO, _require
from bindingmt._private.serialization import json_copy
from bindingmt._private.site_records import (
    capability,
    membership_payload,
    residue_records,
    scope_record,
    sequence,
    source_record,
    text,
    validate,
)
from bindingmt._private.smonitor import signal


def _producer(operation):
    from bindingmt import __version__

    return {
        'library': 'bindingmt',
        'bindingmt_version': __version__,
        'operation': operation,
    }


def _record(source, claim, criterion, name, structures, states, extensions, views):
    return {
        'format': 'bindingmt.binding_site',
        'format_version': '0.1',
        'experimental': True,
        'record_id': str(uuid4()),
        'name': name,
        'context': {'source': source},
        'views': views,
        'definition': {
            'key': 'definition',
            'claim': claim,
            'criterion': criterion,
            'participants': [source['id']],
            'views': [v['key'] for v in views],
            'scope': scope_record(source['id'], structures, states),
            'evidence_refs': [],
            'binding_support': 'unknown',
        },
        'evidence': [],
        'realization_coverage': {'state': 'not_requested'},
        'derivation': {'parents': [], 'kind': None, 'reason': None},
        'extensions': {} if extensions is None else extensions,
        'producer': _producer('declaration'),
    }


@dataclass(frozen=True, init=False, repr=False)
class BindingSite(_InterfaceRecordIO):
    """A declared hypothesis; source addresses are unresolved molecular membership."""

    _json: str

    @arg_digest(config='bindingmt._argdigest')
    def __init__(self, record, *, experimental=False):
        _require(
            experimental is True,
            'experimental',
            'Opt in to experimental site 0.1/0.2/0.3/0.4.',
        )
        snapshot = json_copy(record, 'record')
        validate(snapshot)
        object.__setattr__(self, '_json', json.dumps(snapshot, allow_nan=False))

    @classmethod
    @signal(tags=['api', 'site', 'declaration'])
    @arg_digest(config='bindingmt._argdigest')
    def from_residues(
        cls,
        *,
        source,
        address_space,
        residues,
        name=None,
        role='annotated_members',
        structure_indices=None,
        chemical_state_indices=None,
        extensions=None,
    ):
        """Declare ordered source addresses without coordinates, mapping or evidence."""
        source = source_record(source)
        _require(
            text(address_space) and text(role),
            'address_space',
            'Declare address space and role strings.',
        )
        capability(address_space == 'author_residue', 'site residue address space')
        capability(role == 'annotated_members', 'site annotation role')
        views = [
            {
                'key': 'membership',
                'kind': 'residue_annotation',
                'role': role,
                'scope': scope_record(
                    source['id'], structure_indices, chemical_state_indices
                ),
                'payload': {
                    'address_space': address_space,
                    'entries': residue_records(residues),
                    'resolution': {
                        'state': 'unresolved',
                        'cause': 'no_correspondence_supplied',
                    },
                },
            }
        ]
        record = _record(
            source,
            'Candidate association region; binding support unknown',
            {'kind': 'annotation_membership', 'rule': 'exact_listed_entries'},
            name,
            structure_indices,
            chemical_state_indices,
            extensions,
            views,
        )
        record = json_copy(record, 'record')
        record['content_sha256'] = _digest(record)
        return cls(record, experimental=True)

    @classmethod
    @signal(tags=['api', 'site', 'declaration'])
    @arg_digest(config='bindingmt._argdigest')
    def from_definition(
        cls,
        *,
        source,
        claim,
        criterion,
        name=None,
        structure_indices=None,
        chemical_state_indices=None,
        extensions=None,
    ):
        """Declare an unevaluated criterion, without an invented empty membership."""
        record = _record(
            source_record(source),
            claim,
            criterion,
            name,
            structure_indices,
            chemical_state_indices,
            extensions,
            [],
        )
        record = json_copy(record, 'definition')
        record['content_sha256'] = _digest(record)
        return cls(record, experimental=True)

    @classmethod
    def _from_members(
        cls,
        source,
        indices,
        element,
        name,
        role,
        structures,
        states,
        extensions,
        *,
        captured_selection=None,
    ):
        source = source_record(source)
        _require(text(role), 'role', 'Declare a membership role string.')
        selected = captured_selection is not None
        capability(
            role == ('selected_members' if selected else 'declared_members'),
            'site membership role',
        )
        payload = membership_payload(indices, element)
        if selected:
            payload.update(
                resolution={'state': 'resolved', 'cause': 'captured_selection'},
                selection=captured_selection,
            )
        views = [
            {
                'key': 'membership',
                'kind': 'molecular_membership',
                'role': role,
                'scope': scope_record(source['id'], structures, states),
                'payload': payload,
            }
        ]
        record = _record(
            source,
            'Candidate association region; binding support unknown',
            {
                'kind': 'evaluated_source_selection'
                if selected
                else 'explicit_source_membership',
                'rule': 'fixed_evaluated_membership'
                if selected
                else 'exact_listed_indices',
                'element': element,
            },
            name,
            structures,
            states,
            extensions,
            views,
        )
        record['format_version'] = '0.4' if selected else '0.2'
        if selected:
            record['producer'] = _producer('from_selection')
        record = json_copy(record, 'record')
        record['content_sha256'] = _digest(record)
        return cls(record, experimental=True)

    @classmethod
    @signal(tags=['api', 'site', 'declaration'])
    @arg_digest(config='bindingmt._argdigest')
    def from_atoms(
        cls,
        *,
        source,
        atom_indices,
        name=None,
        role='declared_members',
        structure_indices=None,
        chemical_state_indices=None,
        extensions=None,
    ):
        """Declare source atom indices; no source bound or geometry is checked."""
        return cls._from_members(
            source,
            atom_indices,
            'atom',
            name,
            role,
            structure_indices,
            chemical_state_indices,
            extensions,
        )

    @classmethod
    @signal(tags=['api', 'site', 'declaration'])
    @arg_digest(config='bindingmt._argdigest')
    def from_groups(
        cls,
        *,
        source,
        group_indices,
        name=None,
        role='declared_members',
        structure_indices=None,
        chemical_state_indices=None,
        extensions=None,
    ):
        """Declare source group indices; no residue-ID parsing or atom expansion."""
        return cls._from_members(
            source,
            group_indices,
            'group',
            name,
            role,
            structure_indices,
            chemical_state_indices,
            extensions,
        )

    @classmethod
    @signal(tags=['api', 'site', 'selection'])
    @arg_digest(config='bindingmt._argdigest')
    def from_selection(
        cls,
        molecular_system,
        selection='all',
        *,
        source,
        element='atom',
        syntax='MolSysMT',
        structure_indices=None,
        chemical_state=None,
        name=None,
        role='selected_members',
        extensions=None,
        max_source_atoms=100000,
        max_structure_occurrences=4096,
        max_dense_pairs=1000000,
    ):
        """Capture fixed atom/group membership selected by public MolSysMT.

        Without a structure, query the resident topology/chemical inventory;
        coordinate-requiring expressions need one explicit source structure.
        None state means reference for inventory, structure association otherwise.
        Source association and consistent inputs are caller preconditions. The
        experimental 0.4 Site preserves evaluated context, not a live rule or
        a binding finding. See ``devguide/binding_site_selection.md``.
        """
        from bindingmt._private.site_selection import capture, request_record

        source = source_record(source)
        capability(role == 'selected_members', 'selected site membership role')
        _require(name is None or text(name), 'name', 'Use a nonempty name or None.')
        if extensions is not None:
            extensions = json_copy(extensions, 'extensions')
        request = request_record(
            selection,
            element,
            syntax,
            structure_indices,
            chemical_state,
            max_source_atoms=max_source_atoms,
            max_structure_occurrences=max_structure_occurrences,
            max_dense_pairs=max_dense_pairs,
        )
        indices, captured = capture(molecular_system, request)
        structures = request['structure_indices']
        return cls._from_members(
            source,
            indices,
            element,
            name,
            role,
            structures,
            None
            if structures is None
            else [captured['effective']['chemical_state_index']],
            extensions,
            captured_selection=captured,
        )

    @property
    def site_id(self):
        return self.record_id

    @classmethod
    def _from_geometry(
        cls,
        source,
        kind,
        center,
        frame,
        name,
        role,
        structures,
        states,
        extensions,
        *,
        radius=None,
        lengths=None,
        axes=None,
    ):
        from bindingmt._private.site_geometry import geometry_payload

        source = source_record(source)
        _require(text(role), 'role', 'Declare a geometry role string.')
        capability(role == 'declared_geometry', 'site geometry role')
        views = [
            {
                'key': 'geometry',
                'kind': kind,
                'role': role,
                'scope': scope_record(source['id'], structures, states),
                'payload': geometry_payload(
                    kind, center, frame, radius=radius, lengths=lengths, axes=axes
                ),
            }
        ]
        record = _record(
            source,
            'Candidate association region; binding support unknown',
            {'kind': 'explicit_geometric_support', 'shape': kind},
            name,
            structures,
            states,
            extensions,
            views,
        )
        record['format_version'] = '0.3'
        record = json_copy(record, 'record')
        record['content_sha256'] = _digest(record)
        return cls(record, experimental=True)

    @classmethod
    @signal(tags=['api', 'site', 'geometry'])
    @arg_digest(config='bindingmt._argdigest')
    def from_sphere(
        cls,
        *,
        source,
        center,
        radius,
        frame,
        name=None,
        role='declared_geometry',
        structure_indices=None,
        chemical_state_indices=None,
        extensions=None,
    ):
        """Declare a positive-radius sphere in a stated single-image Cartesian frame."""
        return cls._from_geometry(
            source,
            'sphere',
            center,
            frame,
            name,
            role,
            structure_indices,
            chemical_state_indices,
            extensions,
            radius=radius,
        )

    @classmethod
    @signal(tags=['api', 'site', 'geometry'])
    @arg_digest(config='bindingmt._argdigest')
    def from_box(
        cls,
        *,
        source,
        center,
        lengths,
        frame,
        axes=None,
        name=None,
        role='declared_geometry',
        structure_indices=None,
        chemical_state_indices=None,
        extensions=None,
    ):
        """Declare full positive edge lengths and preserve optional oriented local axes."""
        return cls._from_geometry(
            source,
            'box',
            center,
            frame,
            name,
            role,
            structure_indices,
            chemical_state_indices,
            extensions,
            lengths=lengths,
            axes=axes,
        )

    @property
    def geometry(self):
        """Detached physical quantities/axes, without molecular lookup or placement."""
        from bindingmt._private.site_geometry import public_geometry

        record = self.to_dict()
        if not record['views'] or record['views'][0]['kind'] not in ('sphere', 'box'):
            return None
        view = record['views'][0]
        return public_geometry(view['kind'], view['payload'])

    @property
    def context(self):
        return self.to_dict()['context']

    @property
    def views(self):
        return self.to_dict()['views']

    @property
    def definition(self):
        return self.to_dict()['definition']

    @property
    def evidence(self):
        return self.to_dict()['evidence']

    @property
    def definition_ref(self):
        record = self.to_dict()
        return {
            'record_id': record['record_id'],
            'content_sha256': record['content_sha256'],
            'definition_key': 'definition',
        }

    @signal(tags=['api', 'site', 'projection'])
    @arg_digest(config='bindingmt._argdigest')
    def resolve_region(
        self,
        molecular_system,
        *,
        target=None,
        atom_map=None,
        group_correspondence=None,
        structure_correspondence=None,
        view=None,
        structure_index=None,
        target_structure_index=None,
        chemical_state=None,
        expansion=None,
        coverage_policy='require_complete',
        context_policy='require_match',
        max_source_atoms=100000,
        max_structure_occurrences=4096,
    ):
        """Project source atom/group entries through public MolSysMT queries.

        ``molecular_system`` is the original resident MolSys associated with the
        declaration by the caller. A separate ``target`` requires explicit atom
        and structure maps. Groups require ``source_members`` or
        ``target_group_members`` expansion. Default coverage refuses omissions.
        The detached experimental 0.1 result records checked correspondence;
        it does not assert source authenticity, geometric inclusion or binding.
        See ``devguide/binding_site_region.md`` for scope and trust preconditions.
        """
        from bindingmt._private.site_region import resolve

        return resolve(
            self,
            molecular_system,
            target=target,
            atom_map=atom_map,
            group_correspondence=group_correspondence,
            structure_correspondence=structure_correspondence,
            view=view,
            structure_index=structure_index,
            target_structure_index=target_structure_index,
            chemical_state=chemical_state,
            expansion=expansion,
            coverage_policy=coverage_policy,
            context_policy=context_policy,
            max_source_atoms=max_source_atoms,
            max_structure_occurrences=max_structure_occurrences,
        )

    @staticmethod
    @signal(tags=['api', 'site', 'persistence'])
    @arg_digest(config='bindingmt._argdigest')
    def region_from_dict(record, *, experimental=False):
        """Admit a saved region result without querying a molecular source.

        This checks historical record integrity and correspondence, not current
        applicability. Explicit experimental opt-in is required.
        """
        from bindingmt._private.site_region import read

        return read(record, experimental)

    @signal(tags=['api', 'site', 'application'])
    @arg_digest(config='bindingmt._argdigest')
    def to_search_domain(
        self,
        molecular_system,
        *,
        method,
        view=None,
        structure_index=None,
        chemical_state=None,
        expansion=None,
        padding=None,
        assurance_scope='placed_geometry',
        context_policy='require_match',
        max_source_atoms=100000,
        max_structure_occurrences=4096,
        frame_id=None,
    ):
        """Enclose fixed native membership and retain actual-input use conditions.

        Supply ``method='enclosing_box'`` and explicit length padding. Original
        source association is a caller precondition. The resulting BoxRegion
        retains its Site/projection, geometry and exact position witness. Its
        use requires the actual DockingProblem and checks source and normalized
        receptor positions. It asserts no current binding finding.
        Declared spheres/axis-aligned boxes also admit their named geometry
        method and explicit ``frame_id`` associating the saved frame with the
        original source. This caller association is not an inferred alignment.
        See ``devguide/binding_site_application.md`` and
        ``devguide/binding_site_geometry_application.md`` for admitted profiles.
        """
        from bindingmt._private.site_application import derive

        return derive(
            self,
            molecular_system,
            method=method,
            view=view,
            structure_index=structure_index,
            chemical_state=chemical_state,
            expansion=expansion,
            padding=padding,
            assurance_scope=assurance_scope,
            context_policy=context_policy,
            max_source_atoms=max_source_atoms,
            max_structure_occurrences=max_structure_occurrences,
            frame_id=frame_id,
        )

    def _derive(self, record, kind, reason):
        _require(text(reason), 'reason', 'Record why this declaration is derived.')
        record['record_id'] = str(uuid4())
        record['derivation'] = {
            'parents': [
                {k: v for k, v in self.definition_ref.items() if k != 'definition_key'}
            ],
            'kind': kind,
            'reason': reason,
        }
        record['producer'] = _producer(kind)
        record = json_copy(record, 'record')
        record['content_sha256'] = _digest(record)
        return type(self)(record, experimental=True)

    @signal(tags=['api', 'site', 'interpretation'])
    @arg_digest(config='bindingmt._argdigest')
    def with_definition(self, *, claim, criterion, reason, evidence_refs=()):
        """Revise interpretation; preserve observations and their original claim refs."""
        record = self.to_dict()
        record['definition'].update(
            claim=claim,
            criterion=criterion,
            evidence_refs=sequence(evidence_refs, 'evidence_refs'),
        )
        return self._derive(record, 'definition_revision', reason)

    @signal(tags=['api', 'site', 'annotation'])
    @arg_digest(config='bindingmt._argdigest')
    def with_evidence(
        self,
        *,
        key,
        source,
        observation,
        attribution,
        structure_indices=None,
        chemical_state_indices=None,
    ):
        """Add an attributed annotation and explicitly cite it, without new measurement."""
        record = self.to_dict()
        record['evidence'].append(
            {
                'key': key,
                'kind': 'source_annotation',
                'source': source_record(source),
                'observation': observation,
                'attribution': attribution,
                'participants': record['definition']['participants'],
                'views': record['definition']['views'],
                'scope': scope_record(
                    record['context']['source']['id'],
                    structure_indices,
                    chemical_state_indices,
                ),
                'original_claim_ref': self.definition_ref,
            }
        )
        record['definition']['evidence_refs'].append(key)
        return self._derive(
            record,
            'evidence_addition',
            'Explicitly cite an attributed source annotation',
        )

    def inspect(self):
        """Report declared meaning and missing prerequisites without resolving sources."""
        record = self.to_dict()
        payload = record['views'][0]['payload'] if record['views'] else None
        molecular = (
            payload is not None and record['views'][0]['kind'] == 'molecular_membership'
        )
        geometric = payload is not None and record['views'][0]['kind'] in (
            'sphere',
            'box',
        )
        entries = (
            None
            if payload is None or geometric
            else payload['source_indices']
            if molecular
            else payload['entries']
        )
        membership = {
            'state': 'not_supplied'
            if entries is None
            else 'declared_empty'
            if not entries
            else 'declared'
            if molecular
            else 'unresolved',
            'entry_count': None if entries is None else len(entries),
        }
        if molecular:
            membership.update(
                element=payload['element'], index_space=payload['index_space']
            )
        report = {
            'record_id': self.record_id,
            'name': record['name'],
            'definition_ref': self.definition_ref,
            'source': record['context']['source'],
            'definition_scope': record['definition']['scope'],
            'membership': membership,
            'binding_support': 'unknown',
            'evidence_count': len(record['evidence']),
            'current_evidence_refs': record['definition']['evidence_refs'],
            'evidence_scopes': {e['key']: e['scope'] for e in record['evidence']},
            'realization_coverage': record['realization_coverage'],
            'molecular_resolution': 'not_performed',
            'asset_access': 'not_checked',
            'capabilities': {
                'declaration': True,
                'persistence': True,
                'molecular_projection': molecular,
                'domain_placement': molecular,
            },
        }
        if geometric:
            report['capabilities']['domain_placement'] = record['views'][0][
                'kind'
            ] == 'sphere' or self.geometry['axes'].tolist() == [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
            report['geometry'] = {
                'state': 'declared',
                'kind': record['views'][0]['kind'],
                'frame': payload['frame'],
                'placement_assurance': 'not_assessed',
            }
        if molecular and 'selection' in payload:
            report['selection_capture'] = payload['selection']
        return report

    def describe(self):
        report = self.inspect()
        representation = (
            f'{report["geometry"]["kind"]} geometry in frame {report["geometry"]["frame"]["id"]}'
            if 'geometry' in report
            else f'{report["membership"]["state"]} membership'
        )
        return (
            f'BindingSite {report["name"] or report["record_id"]}: {representation}; '
            f'{report["evidence_count"]} attributed observations. '
            'Binding support is unknown; molecular resolution is not performed.'
        )
