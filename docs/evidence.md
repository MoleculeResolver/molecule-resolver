# Evidence Output (`include_evidence`)

`find_single_molecule_crosschecked(..., include_evidence=True)` returns a
`ResolutionResult` instead of only a `Molecule`.

## Response fields

- `best_molecule`: Selected molecule (or `None` for unresolved ties).
- `ranked_candidates`: List of `CandidateEvidence` entries sorted by score.
- `grouped_by_structure`: Raw grouped molecules by SMILES.
- `selected_smiles`: SMILES selected as best candidate.
- `selection_reason`: Why this result form was returned.

## CandidateEvidence fields

- `service_agreement_count`: Number of unique supporting services.
- `identifier_concordance_count`: Number of unique supporting identifiers.
- `synonym_overlap_count`: Repeated synonym matches across supporting molecules.
- `score_breakdown`: Explicit points per evidence component.
- `total_score`: Sum of `score_breakdown`.
