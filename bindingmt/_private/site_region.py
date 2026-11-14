"""Bounded native site projection, with detached historical result admission.

MolSysMT owns molecular queries and typed atom maps. This module arranges the
site's declared entries and checks their recorded correspondence, not chemistry,
asset identity, geometric inclusion or molecular reconstruction.
"""

import json
from numbers import Integral
from uuid import UUID, uuid4

from bindingmt._private.interface_records import _digest, _require
from bindingmt._private.serialization import json_copy
from bindingmt._private.site_records import keys
from bindingmt._private.smonitor import CapabilityMismatchError

_REQUEST = {
    'derived',
    'view',
    'structure_index',
    'target_structure_index',
    'chemical_state',
    'expansion',
    'coverage_policy',
    'context_policy',
    'max_source_atoms',
    'max_structure_occurrences',
}
_COUNTS = ('n_atoms', 'n_groups', 'n_structures', 'n_chemical_states')
_OBSERVATION = {*_COUNTS, 'atom_group_indices', 'chemical_state_index', 'group_members'}
_TRUST = {
    'caller_preconditions': [
        'original source asset and revision association',
        'consistent source and target throughout check and use',
        'derived state collection preserves source positions and meaning',
    ],
    'checked': [
        'counts',
        'source selection',
        'atom parents',
        'selected state positions',
        'group expansion',
    ],
    'not_assessed': [
        'source authenticity',
        'coordinate frame',
        'geometric inclusion',
        'chemical equivalence',
        'current binding finding',
        'binding support',
        'consumer execution',
    ],
}


def _supported(condition, name):
    if not condition:
        raise CapabilityMismatchError(capability=name, engine='bindingmt.site.region')


def _same(left, right):
    return json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)


def _index(value, field, maximum=None, minimum=0):
    _require(
        isinstance(value, Integral)
        and not isinstance(value, bool)
        and minimum <= value
        and (maximum is None or value < maximum),
        field,
        'Use an integer in the admitted axis bounds.',
    )
    return int(value)


def _vector(value, field, length, maximum, minimum=0):
    _require(
        isinstance(value, list) and len(value) == length,
        field,
        'Preserve the complete ordered correspondence axis.',
    )
    return [_index(i, field, maximum, minimum) for i in value]


def _request(request):
    keys(request, _REQUEST, 'request')
    _require(
        type(request['derived']) is bool, 'target', 'Record an explicit target mode.'
    )
    for field in ('max_source_atoms', 'max_structure_occurrences'):
        _index(request[field], field, minimum=1)
    for field in ('structure_index', 'target_structure_index'):
        if request[field] is not None:
            _index(request[field], field)
    state = request['chemical_state']
    if state is not None and not (isinstance(state, str) and state == 'structure'):
        _index(state, 'chemical_state')
    _supported(request['context_policy'] == 'require_match', 'region context policy')
    _supported(
        request['coverage_policy'] in ('require_complete', 'allow_partial'),
        'region coverage policy',
    )


def _view(site, request):
    views = [
        v
        for v in site.views
        if v['kind'] == 'molecular_membership'
        and v['payload']['element'] in ('atom', 'group')
    ]
    if request['view'] is not None:
        views = [v for v in views if v['key'] == request['view']]
    _supported(len(views) == 1, 'single native atom/group membership view')
    view = views[0]
    expansion = request['expansion']
    if view['payload']['element'] == 'atom':
        _supported(expansion in (None, 'listed_atoms'), 'atom region expansion')
        expansion = 'listed_atoms'
    else:
        _supported(
            expansion in ('source_members', 'target_group_members'),
            'explicit group region expansion',
        )
    return view, expansion


def _counts(counts, request):
    for name in _COUNTS:
        limit = (
            request['max_source_atoms']
            if name in ('n_atoms', 'n_groups')
            else request['max_structure_occurrences']
        )
        _index(counts[name], name, limit + 1)
    _require(
        counts['n_structures'] > 0 and counts['n_chemical_states'] > 0,
        'molecular_system',
        'Require structures and a known chemical-state collection.',
    )


def _axes(view, request, source, target, structures):
    _counts(source, request)
    _counts(target, request)
    scope = view['scope']
    declared = scope.get('source_structure_indices', [])
    _require(
        len(declared) <= request['max_structure_occurrences'],
        'scope',
        'Declared occurrences exceed the admitted budget.',
    )
    for i in declared:
        _index(i, 'scope', source['n_structures'])
    for state in scope.get('chemical_state_indices', []):
        if state is not None:
            _index(state, 'scope', source['n_chemical_states'])
    chosen = request['structure_index']
    if chosen is None:
        distinct = set(declared)
        _require(
            len(distinct) == 1 or not distinct and source['n_structures'] == 1,
            'structure_index',
            'Choose one original source structure explicitly.',
        )
        chosen = next(iter(distinct)) if distinct else 0
    chosen = _index(chosen, 'structure_index', source['n_structures'])
    _require(
        not declared or chosen in declared,
        'structure_index',
        'The chosen source structure must belong to the declared scope.',
    )
    if request['derived']:
        structures = _vector(
            structures,
            'structure_correspondence',
            target['n_structures'],
            source['n_structures'],
        )
        _require(
            source['n_chemical_states'] == target['n_chemical_states'],
            'chemical_state',
            'Derived state collection cardinality must be preserved.',
        )
    else:
        _require(
            structures is None and source == target,
            'target',
            'Direct mode uses exactly the same source axis.',
        )
        structures = list(range(source['n_structures']))
    inverse = [[] for _ in range(source['n_structures'])]
    for i, origin in enumerate(structures):
        inverse[origin].append(i)
    target_index = request['target_structure_index']
    if target_index is None:
        _require(
            len(inverse[chosen]) == 1,
            'target_structure_index',
            'Choose one target occurrence; source references can repeat or be absent.',
        )
        target_index = inverse[chosen][0]
    target_index = _index(
        target_index, 'target_structure_index', target['n_structures']
    )
    _require(
        structures[target_index] == chosen,
        'target_structure_index',
        'Target occurrence must refer to the chosen source structure.',
    )
    return {
        'source_structure_index': chosen,
        'target_structure_index': target_index,
        'target_to_source': structures,
        'source_to_target': inverse,
        'declaration_occurrences': [
            i for i, origin in enumerate(declared) if origin == chosen
        ],
    }


def _portable(value):
    from bindingmt.core.search_domain import _descriptor

    return _descriptor(value)


def _map(raw, source, target, request):
    import molsysmt as msm

    if not request['derived']:
        _require(raw is None, 'atom_map', 'Direct source mode needs no external map.')
        raw = dict(
            schema='molsysmt.atom_map@1',
            n_source_atoms=source['n_atoms'],
            n_output_atoms=source['n_atoms'],
            source_to_output_atom_indices=list(range(source['n_atoms'])),
            output_to_source_atom_indices=list(range(source['n_atoms'])),
        )
    _require(
        isinstance(raw, dict),
        'atom_map',
        'An explicit owner atom-map report is required.',
    )

    # Bound stored stages and total map data before owner axis allocation.
    data_items = 0

    def budget(value, depth=0):
        nonlocal data_items
        if isinstance(value, (dict, list)):
            data_items += len(value)
        _require(
            data_items
            <= 32
            * (request['max_source_atoms'] + request['max_structure_occurrences']),
            'atom_map',
            'Total atom-map data exceeds the admitted budget.',
        )
        _require(
            depth <= 32, 'atom_map', 'Atom-map nesting exceeds the admitted budget.'
        )
        if isinstance(value, dict):
            for name, item in value.items():
                if name in ('n_source_atoms', 'n_output_atoms'):
                    _index(item, 'atom_map', request['max_source_atoms'] + 1)
                budget(item, depth + 1)
        elif isinstance(value, list):
            _require(
                len(value)
                <= max(
                    request['max_source_atoms'], request['max_structure_occurrences']
                ),
                'atom_map',
                'Atom-map arrays exceed the admitted budget.',
            )
            for item in value:
                if isinstance(item, (dict, list)):
                    budget(item, depth + 1)

    budget(raw)
    canonical = _portable(msm.compose_atom_maps([raw]))
    _require(
        canonical['n_source_atoms'] == source['n_atoms']
        and canonical['n_output_atoms'] == target['n_atoms'],
        'atom_map',
        'Map counts must match the associated source and target axes.',
    )
    _supported(
        canonical['source_to_representative_atom_indices']
        == canonical['source_to_output_atom_indices'],
        'unmerged original atom identity',
    )
    return canonical


def _parents(observation):
    return _vector(
        observation['atom_group_indices'],
        'atom_group_indices',
        observation['n_atoms'],
        observation['n_groups'],
        -1,
    )


def _groups(source, target, mapping, declared, derived, element):
    left, right = _parents(source), _parents(target)
    if element == 'group':
        _supported(-1 not in left and -1 not in right, 'known group parentage')
    witnesses, reverse = {}, {}
    for atom, output in enumerate(mapping['source_to_output_atom_indices']):
        if output < 0 or left[atom] < 0 or right[output] < 0:
            continue
        a, b = left[atom], right[output]
        _supported(a not in witnesses or witnesses[a] == b, 'unsplit group identity')
        _supported(b not in reverse or reverse[b] == a, 'unmerged group identity')
        witnesses[a], reverse[b] = b, a
    if not derived:
        _require(
            declared is None,
            'group_correspondence',
            'Direct source mode needs no external group map.',
        )
        return list(range(source['n_groups']))
    if declared is None:
        return [witnesses.get(i) for i in range(source['n_groups'])]
    result = _vector(
        declared, 'group_correspondence', source['n_groups'], target['n_groups'], -1
    )
    present = [i for i in result if i >= 0]
    _require(
        len(set(present)) == len(present)
        and all(result[a] == b for a, b in witnesses.items()),
        'group_correspondence',
        'Explicit injective group map must agree with every surviving atom witness.',
    )
    return result


def _layout(site, request, observations, correspondence):
    _request(request)
    view, expansion = _view(site, request)
    keys(observations, {'source', 'target'}, 'observations')
    source, target = observations['source'], observations['target']
    for item in (source, target):
        keys(item, _OBSERVATION, 'observations')
        _counts(item, request)
    _require(
        request['derived'] or _same(source, target),
        'observations',
        'Direct mode must preserve the same source observations on both axes.',
    )
    keys(
        correspondence,
        {
            'atom_map',
            'canonical_atom_map',
            'group_correspondence',
            'structure_correspondence',
        },
        'correspondence',
    )

    def counts(obs):
        return {k: obs[k] for k in _COUNTS}

    axes = _axes(
        view,
        request,
        counts(source),
        counts(target),
        correspondence['structure_correspondence'],
    )
    state = _index(
        source['chemical_state_index'], 'chemical_state', source['n_chemical_states']
    )
    _require(
        _index(
            target['chemical_state_index'],
            'chemical_state',
            target['n_chemical_states'],
        )
        == state,
        'chemical_state',
        'Selected source and target state positions must agree.',
    )
    requested_state = request['chemical_state']
    if requested_state is not None and requested_state != 'structure':
        _require(
            requested_state == state,
            'chemical_state',
            'Explicit state must match the selected structure association.',
        )
    for occurrence in axes['declaration_occurrences']:
        declared_state = view['scope']['chemical_state_indices'][occurrence]
        _require(
            declared_state is None or declared_state == state,
            'chemical_state',
            'Selected known declaration states must match the actual association.',
        )
    mapping = _map(correspondence['atom_map'], source, target, request)
    if correspondence['canonical_atom_map'] is not None:
        # Software version is historical producer information, not map semantics.
        stored = _map(
            correspondence['canonical_atom_map'],
            source,
            target,
            {**request, 'derived': True},
        )
        _require(
            {k: v for k, v in stored.items() if k != 'software'}
            == {k: v for k, v in mapping.items() if k != 'software'},
            'atom_map',
            'Stored canonical map must agree with the original owner report.',
        )
    element = view['payload']['element']
    groups = _groups(
        source,
        target,
        mapping,
        correspondence['group_correspondence'],
        request['derived'],
        element,
    )
    indices = view['payload']['source_indices']
    _require(
        len(indices) <= request['max_source_atoms'],
        'membership',
        'Membership exceeds the admitted budget.',
    )
    for i in indices:
        _index(i, 'membership', source['n_atoms' if element == 'atom' else 'n_groups'])
    if element == 'group':
        _require(
            all(groups[i] is not None for i in indices),
            'group_correspondence',
            'Groups without surviving atom witnesses require explicit correspondence, including empty groups.',
        )
    axes['chemical_state_index'] = state
    return view, expansion, axes, mapping, groups


def _projection(view, expansion, axes, mapping, groups, observations, request):
    source, target = observations['source'], observations['target']
    element, indices = view['payload']['element'], view['payload']['source_indices']
    forward, inverse = (
        mapping['source_to_output_atom_indices'],
        mapping['output_to_source_atom_indices'],
    )
    for obs in (source, target):
        _require(
            isinstance(obs['group_members'], list)
            and len(obs['group_members'])
            == (len(indices) if element == 'group' else 0),
            'group_members',
            'Keep one expansion observation per declared group entry.',
        )
    # Linear historical consistency checks; native expansion is owner supplied.
    expected_source, expected_target = {}, {}
    if element == 'group':
        expected_source = {i: [] for i in indices}
        expected_target = {groups[i]: [] for i in indices if groups[i] >= 0}
        for atom, parent in enumerate(source['atom_group_indices']):
            if parent in expected_source:
                expected_source[parent].append(atom)
        for atom, parent in enumerate(target['atom_group_indices']):
            if parent in expected_target:
                expected_target[parent].append(atom)
    entries, exported = [], []
    for position, index in enumerate(indices):
        group = groups[index] if element == 'group' else None
        originals = source['group_members'][position] if element == 'group' else [index]
        if element == 'group':
            expected = expected_source[index]
            _require(
                _same(originals, expected),
                'group_members',
                'Source expansion must agree with the complete recorded parent axis.',
            )
            actual = target['group_members'][position]
            expected = None if group < 0 else expected_target[group]
            _require(
                _same(actual, expected),
                'group_members',
                'Target expansion must agree with the complete recorded parent axis.',
            )
        missing = [i for i in originals if forward[i] < 0]
        omitted_group = element == 'group' and group < 0
        _require(
            request['coverage_policy'] == 'allow_partial'
            or not missing
            and not omitted_group,
            'coverage_policy',
            'Complete coverage refuses known omitted original members.',
        )
        output = (
            [forward[i] for i in originals if forward[i] >= 0]
            if expansion != 'target_group_members'
            else []
            if group < 0
            else actual
        )
        generated = [i for i in output if inverse[i] < 0]
        entries.append(
            {
                'source_index': index,
                'target_group_index': group,
                'source_atom_indices': originals,
                'target_atom_indices': output,
                'output_source_atom_indices': [
                    None if inverse[i] < 0 else inverse[i] for i in output
                ],
                'missing_source_atom_indices': missing,
                'generated_target_atom_indices': generated,
                'omitted_group': omitted_group,
                'state': 'omitted'
                if omitted_group
                else 'partial'
                if missing
                else 'empty'
                if not output
                else 'complete',
            }
        )
        exported.extend(output)
    return {
        'view_key': view['key'],
        'element': element,
        'expansion': expansion,
        'scope': view['scope'],
        'axes': axes,
        'group_correspondence': groups,
        'entries': entries,
        'exported_target_atom_indices': list(dict.fromkeys(exported)),
        'coverage': 'declared_empty'
        if not entries
        else 'partial'
        if any(e['missing_source_atom_indices'] or e['omitted_group'] for e in entries)
        else 'complete',
        'assurance': _TRUST,
    }


def _counts_from_source(msm, molecular_system, request):
    _supported(
        msm.get_form(molecular_system) == 'molsysmt.MolSys', 'resident native MolSys'
    )
    counts = {
        name: _index(msm.get(molecular_system, **{name: True}), name)
        for name in _COUNTS
    }
    _counts(counts, request)
    return counts


def _observe(msm, system, counts, structure):
    parents = msm.get(
        system,
        element='atom',
        group_index=True,
        structure_indices=[structure],
        chemical_state='structure',
    )
    # MolSysMT returns nullable parentage; -1 is only a stored unknown sentinel.
    parents = [
        int(i) if isinstance(i, Integral) and not isinstance(i, bool) else -1
        for i in parents
    ]
    states = msm.get(
        system, structure_indices=[structure], structure_chemical_state_index=True
    )
    state = states[0]
    _require(
        isinstance(state, Integral) and not isinstance(state, bool),
        'chemical_state',
        'Require a known selected structure-state association.',
    )
    return {
        **counts,
        'atom_group_indices': parents,
        'chemical_state_index': int(state),
        'group_members': [],
    }


def _expand(msm, system, groups, structure):
    # A single ordered owner query; never send an omission sentinel as an index.
    selected = [i for i in groups if i >= 0]
    rows = (
        msm.get(
            system,
            element='group',
            selection=selected,
            atom_index=True,
            structure_indices=[structure],
            chemical_state='structure',
        )
        if selected
        else []
    )
    rows = iter(_portable(rows))
    return [next(rows) if i >= 0 else None for i in groups]


def resolve(
    site,
    molecular_system,
    *,
    target,
    atom_map,
    group_correspondence,
    structure_correspondence,
    **request,
):
    import molsysmt as msm

    from bindingmt.core.site import _producer

    request = _portable({'derived': target is not None, **request})
    _request(request)
    view, expansion = _view(site, request)
    target = molecular_system if target is None else target
    # Do not copy an entire molecule, or infer source identity from counts/IDs.
    source_counts = _counts_from_source(msm, molecular_system, request)
    target_counts = (
        source_counts
        if not request['derived']
        else _counts_from_source(msm, target, request)
    )
    structures = _portable(structure_correspondence)
    axes = _axes(view, request, source_counts, target_counts, structures)
    observations = {
        'source': _observe(
            msm, molecular_system, source_counts, axes['source_structure_index']
        ),
        'target': _observe(msm, target, target_counts, axes['target_structure_index']),
    }
    correspondence = {
        'atom_map': _portable(atom_map),
        'canonical_atom_map': None,
        'group_correspondence': _portable(group_correspondence),
        'structure_correspondence': structures,
    }
    view, expansion, axes, mapping, groups = _layout(
        site, request, observations, correspondence
    )
    indices = view['payload']['source_indices']

    def selected_members():
        return _portable(
            msm.select(
                molecular_system,
                element=view['payload']['element'],
                selection=indices,
                structure_indices=[axes['source_structure_index']],
                chemical_state='structure',
            )
        )

    _require(
        selected_members() == indices,
        'membership',
        'Owner selection must preserve the declared entry order.',
    )
    if view['payload']['element'] == 'group':
        observations['source']['group_members'] = _expand(
            msm, molecular_system, indices, axes['source_structure_index']
        )
        missing_atoms = any(
            mapping['source_to_output_atom_indices'][i] < 0
            for row in observations['source']['group_members']
            for i in row
        )
        _require(
            request['coverage_policy'] == 'allow_partial'
            or not missing_atoms
            and all(groups[i] >= 0 for i in indices),
            'coverage_policy',
            'Complete coverage refuses omissions before indexed target expansion.',
        )
        observations['target']['group_members'] = _expand(
            msm, target, [groups[i] for i in indices], axes['target_structure_index']
        )
    projection = _projection(
        view, expansion, axes, mapping, groups, observations, request
    )
    # Repeat the same bounded owner observations; no coordinate/snapshot claim.
    for label, system, index in (
        ('source', molecular_system, axes['source_structure_index']),
        ('target', target, axes['target_structure_index']),
    ):
        again = _observe(msm, system, _counts_from_source(msm, system, request), index)
        if view['payload']['element'] == 'group':
            selected = view['payload']['source_indices']
            if label == 'target':
                selected = [groups[i] for i in selected]
            again['group_members'] = _expand(msm, system, selected, index)
        _require(
            again == observations[label],
            'molecular_system',
            'Observed source/target witnesses changed during projection.',
        )
    _require(
        selected_members() == indices,
        'membership',
        'Owner selection changed during projection.',
    )
    correspondence['canonical_atom_map'] = mapping
    record = {
        'format': 'bindingmt.site_region_projection',
        'format_version': '0.1',
        'experimental': True,
        'record_id': str(uuid4()),
        'site': site.to_dict(),
        'definition_ref': site.definition_ref,
        'view_ref': {'key': view['key'], 'site_record_id': site.record_id},
        'request': request,
        'observations': observations,
        'correspondence': correspondence,
        'projection': projection,
        'producer': _producer('resolve_region'),
    }
    record = json_copy(record, 'record')
    record['content_sha256'] = _digest(record)
    return record


def read(record, experimental):
    from bindingmt.core.site import BindingSite

    _require(
        experimental is True, 'experimental', 'Opt in to experimental site-region 0.1.'
    )
    record = json_copy(record, 'record')
    keys(
        record,
        {
            'format',
            'format_version',
            'experimental',
            'record_id',
            'site',
            'definition_ref',
            'view_ref',
            'request',
            'observations',
            'correspondence',
            'projection',
            'producer',
            'content_sha256',
        },
        'record',
    )
    _supported(
        record['format'] == 'bindingmt.site_region_projection'
        and record['format_version'] == '0.1',
        'site-region record version',
    )
    _require(
        record['experimental'] is True and record['content_sha256'] == _digest(record),
        'record',
        'Preserve experimental record integrity.',
    )
    try:
        UUID(record['record_id'])
    except (ValueError, TypeError, AttributeError) as exc:
        _require(False, 'record_id', f'Use a UUID ({type(exc).__name__}).')
    keys(record['producer'], {'library', 'bindingmt_version', 'operation'}, 'producer')
    _require(
        record['producer']['library'] == 'bindingmt'
        and record['producer']['operation'] == 'resolve_region'
        and isinstance(record['producer']['bindingmt_version'], str),
        'producer',
        'Preserve the original projection producer.',
    )
    site = BindingSite.from_dict(record['site'], experimental=True)
    view, expansion, axes, mapping, groups = _layout(
        site, record['request'], record['observations'], record['correspondence']
    )
    _require(
        record['correspondence']['canonical_atom_map'] is not None,
        'atom_map',
        'Preserve the checked owner map.',
    )
    expected = _projection(
        view,
        expansion,
        axes,
        mapping,
        groups,
        record['observations'],
        record['request'],
    )
    _require(
        record['definition_ref'] == site.definition_ref
        and record['view_ref'] == {'key': view['key'], 'site_record_id': site.record_id}
        and _same(record['projection'], expected),
        'projection',
        'Preserve exact site references and ordered projection relationships.',
    )
    return record
