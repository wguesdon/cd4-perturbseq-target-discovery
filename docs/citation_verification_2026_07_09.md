# Citation verification, 2026-07-09

Every reference in `reports/references.bib` was verified by fetching the paper and confirming that it
supports the specific claim it is cited for. 24 of 25 candidate citations passed. The verification was
run before any manuscript text was written, so no claim rests on an unchecked reference.

This file records what failed and what must change as a result.

EDITORIAL NOTE

Excluded / flagged citations:

1. manguso2017 — EXISTS but DOES NOT SUPPORT THE CLAIM AS WRITTEN. The reference is bibliographically correct (Manguso et al., Nature 547(7664):413-418, 2017, doi:10.1038/nature23270), but it is about mouse tumour-cell-intrinsic Ptpn2 loss enhancing sensitivity to T-cell-mediated immunotherapy. It contains no human PTPN2 genetics, no T-cell-intrinsic loss-of-function analysis, and no autoimmunity content. It therefore cannot support the claim that "loss of function of PTPN2 is associated with autoimmunity rather than immunodeficiency." What to do: split the compound claim. Keep manguso2017 only for the "PTPN2 is a negative regulator of IFN-gamma signalling / immunotherapy target" clause, and note that even there the effect it shows is tumour-cell-intrinsic, not T-cell-intrinsic. For the autoimmunity clause, cite instead (a) Wiede et al., "T cell protein tyrosine phosphatase attenuates T cell signaling to maintain tolerance in mice," J Clin Invest 2011;121(12):4758-4774, doi:10.1172/JCI59492 (T-cell PTPN2 deletion causes hyperactivation and autoimmunity), plus (b) a human GWAS source for PTPN2 autoimmune associations, e.g., Todd et al., Nat Genet 2007;39:857-864, doi:10.1038/ng2068, or WTCCC 2007, Nature 447:661-678, doi:10.1038/nature05911. If a single T-cell-intrinsic + autoimmunity citation is preferred, Wiede 2011 alone is the closest fit; otherwise drop manguso2017 from that sentence entirely.

Special-attention items (both included above, with caveats the author should carry into the text):

2. arce2025 — INCLUDED. Your data-provenance file flagged this as unconfirmed, but the fetch verification resolves the flag: the paper was confirmed against the PMC full text (PMC11754113) and the Nature record, doi:10.1038/s41586-024-08314-y, and it supports the claim (pooled CRISPR TF/chromatin-modifier screens across resting Teff, resting Treg, and restimulated Teff). Two caveats to reflect in the manuscript: (a) Year ambiguity — published online 18 Dec 2024, print volume 637 dated 2025; both are defensible and the key/entry use 2025 to match Nature's official citation, so do not "correct" it to 2024. (b) The word "fitness" is loose — the primary readout is IL-2Ra surface expression sorted by FACS (an activation-state marker screen), not a proliferation/dropout fitness screen. If the sentence means a dropout screen specifically, reword to "CRISPR marker (IL-2Ra) screens across resting and stimulated Teff and Treg."

3. manguso2017 is the only entry that failed outright; the remaining 24 all verified by fetch and support their claims, so they are retained. None failed for non-existence, and none failed for being unfetchable.

Optional precision upgrades (entries kept, but the manuscript may want a more exact primary source for one specific sub-claim):
- hart2017 correctly anchors the CEGv2 core-essential set, but the NEGv1 nonessential set it uses originates in Hart et al. 2014, Mol Syst Biol 10(7):733, doi:10.15252/msb.20145216 — add that if you cite hart2017 specifically for NEGv1.
- karczewski2020 is correct for LOEUF and for gnomAD constraint metrics as released, but the pLI/pRec framework was first defined in Lek et al. 2016, Nature 536:285-291, doi:10.1038/nature19057 — cite Lek 2016 if a primary source for pRec is needed.
