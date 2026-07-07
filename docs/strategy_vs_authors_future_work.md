# Our strategy vs the authors' stated future work

Purpose: confirm our novelty claims do not collide with what the source paper's authors say
is next, and name the genuine white space. Source quotes and details are in
[`literature/source_paper_future_work.md`](./literature/source_paper_future_work.md).

## The authors' future work, in five items

1. More guides per target, plus better off-target and donor-vs-batch deconvolution.
2. Single-cell distribution-aware analyses beyond pseudobulk.
3. Polarizing cytokine conditions and new contexts and modalities.
4. Scaling to hundreds of millions of cells.
5. Tighter integration with population genetics and cell atlases.

None of these is a drug-target-selection metric or a known-drug benchmark.

## Mapping our strategy against it

| Our element | Paper status | In authors' future work | Verdict |
| --- | --- | --- | --- |
| Therapeutic-window target-selection score | context-specificity shown descriptively only | no | White space. Safe headline. |
| Immunomodulator benchmark, precision and recall | no benchmarking vs known drug targets | no | White space. |
| In-silico safety-liability axis | not done | no | White space. |
| Druggability filter | never emphasized | no | Keep as a filter. |
| Autoimmune genetics annotation | already done via GWAS, OpenTargets, OneK1K, UKB LoF | partially | Build on theirs. Cite. Do not claim as novel. |
| Th1 and Th17 program scoring | non-polarized culture only | yes | State the non-polarized limitation. |

## Bottom line

Our contribution is a validated target-triage decision layer. It sits orthogonal to the
authors' technology-scaling roadmap, so the novelty position is strong. Two reframes follow:

1. The genetics annotation builds on the authors' own GWAS and OpenTargets work. Cite it,
   do not present it as new.
2. The Th1 and Th17 signals are proxies read from marker genes within stimulated cells,
   because the data has no polarizing conditions. State this as a limitation.
