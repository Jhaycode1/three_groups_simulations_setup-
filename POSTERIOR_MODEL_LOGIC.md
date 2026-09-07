# Three-Groups Posterior Model Logic

## Model Structure: Two Binary Indicators

The NIMBLE model uses **two binary indicator variables** to define three mutually exclusive groups:

```
┌─────────────────────────────────────────────────────────────────┐
│                  Gene Classification                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                   │
│  not_null = 0  →  **NULL GROUP**                                 │
│                   (No effect at all, gene is inactive)            │
│                                                                   │
│  not_null = 1  →  Gene is ACTIVE (has some effect)               │
│    ├─ not_ben = 0  →  **BENEFICIAL GROUP**                       │
│    │                  (Positive effect on the phenotype)         │
│    │                                                              │
│    └─ not_ben = 1  →  **DELETERIOUS GROUP**                      │
│                       (Negative effect on the phenotype)         │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

## MCMC Sampling

During MCMC, each gene and each iteration gets a binary assignment:
- **Iteration 1**: Gene 1 gets `not_null=1, not_ben=0` → **BENEFICIAL**
- **Iteration 2**: Gene 1 gets `not_null=1, not_ben=0` → **BENEFICIAL**
- **Iteration 3**: Gene 1 gets `not_null=0, not_ben=?` → **NULL** (not_ben value is ignored)
- ...and so on

## Posterior Probability Calculation

The posterior probabilities are the **proportion of MCMC iterations** where each gene falls into each group:

```
P(Null | data) = Proportion of iterations where not_null = 0
                = 1 - mean(not_null)

P(Beneficial | data) = Proportion of iterations where (not_null = 1 AND not_ben = 0)
                     = mean((1 - not_ben) * not_null)
                     = mean of [1 if active & beneficial, 0 otherwise]

P(Deleterious | data) = Proportion of iterations where (not_null = 1 AND not_ben = 1)
                      = mean(not_ben * not_null)
                      = mean of [1 if active & deleterious, 0 otherwise]
```

## Why `not_null` is Critical

The `not_null` indicator **controls membership in the active groups**:

1. **Acts as a filter**: When `not_null = 0`, the gene is automatically NULL
   - The `not_ben` value becomes irrelevant
   - Both beneficial and deleterious probabilities are suppressed

2. **Enables joint probability**: The multiplication `(1 - not_ben) * not_null` ensures:
   - A gene can ONLY be beneficial if BOTH conditions are true:
     - `not_null = 1` (the gene is active)
     - `not_ben = 0` (and it's beneficial, not deleterious)
   
3. **Creates three mutually exclusive groups**: The three probabilities always sum to 1.0
   - These are the only three possible outcomes
   - A gene cannot be both null and beneficial
   - A gene cannot be both beneficial and deleterious

## Example: Gene Simulation

For a gene with:
- 200 iterations as NULL
- 80 iterations as BENEFICIAL  
- 20 iterations as DELETERIOUS
- Total 300 iterations

```
P(Null | data)        = 200/300 = 0.667
P(Beneficial | data)  = 80/300  = 0.267
P(Deleterious | data) = 20/300  = 0.067
Sum                   = 1.000 ✓
```

## Summary: What's the Effect of `not_null`?

`not_null` **gates the entire "active gene" state**:
- It determines whether a gene participates in either beneficial OR deleterious categories
- Without `not_null = 1`, a gene cannot be beneficial or deleterious—it must be null
- This creates a hierarchical model: first decide if gene is active, then decide its effect type
- It ensures the three groups are exhaustive and mutually exclusive
