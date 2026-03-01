# Consensus Scoring

`find_single_molecule_crosschecked(..., resolution_mode="consensus")` uses an
explicit score model for grouped structure candidates.

## Score Components

- `crosscheck_points`: `crosscheck_count * 100`
- `opsin_bonus`: `30` if at least one supporting OPSIN name match exists
- `stereo_specificity_bonus`: RDKit-derived stereochemistry signal:
  - defined chiral centers: `20` points each
  - stereochemical bonds: `15` points each
  - unresolved chiral centers: `5` points each

## Deterministic Tie-Breaking

When total score ties, the selector orders by:

1. Higher crosscheck count
2. Higher OPSIN bonus
3. Higher stereo specificity bonus
4. Higher defined chirality bonus
5. Higher bond stereo bonus
6. Lexicographical SMILES fallback

This removes raw SMILES length from score weighting while preserving a stable
deterministic result.
