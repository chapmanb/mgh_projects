## Overview

Looking at homeobox regulatory proteins in two species -- Drosophila and
mammals. The active regulatory complex being studied has 4 members:

| **Drosophila** | **Mammals**
|----------------|---------------------
| Ph             | Ph
| dRing          | Ring1a, Ring1b
| Pc             | Cbx2/M33, Cb4,6,7,8
| Psc, Suz2      | BmiI, Mel18

Psc is Drosophila consists of a 5' set of domains and 3' region rich in
charged amino acids. The homolog of Psc, BmiI contains only the 5' region and
not the 3' charged region. BmiI does not have the activity of Psc.

So the search is on for the homolog in the mammalian complex that does the
work of the 3' Psc region. Pc homolog Cbx2/M33 has Psc activity and is also
highly charged. Deletions of Cbx2 show a correlation between protein charge
and activity.

## Scripts
- `interpro_domain_summary.py` -- Download domains from InterPro and UniProt,
  parse files, and prepare local key/value database.

  - If first time or modified organisms, build HMM search information:

          <pre><code>
          % hmmsearch Chromo_ls.hmm IPR000953.fa > IPR000953-domains.hmmsearch
          % hmmsearch -E 10 cbx2_Pc-no_chromo-align.hmm IPR000953.fa > IPR000953-filter1.hmmsearch
          % hmmsearch zf-C3HC4_ls.hmm IPR001841.fa > IPR001841-domains.hmmsearch
          % hmmsearch -E 10 bmi1_psc-no_zr-align.hmm IPR001841.fa > IPR001841-filter1.hmmsearch
          </pre></code>

  - Then re-run `interpro_domain_summary.py`

- `family_value_cluster.py` -- Perform k-means clustering and write out HTML files with
  results for display.

- `extract_clusters_of_interest.py` -- Renames the cluster files based on
  associated members.

- `tax_data_display.py` -- Prepare a summary taxonomy tree of classifications
  from the k-means clustering.


