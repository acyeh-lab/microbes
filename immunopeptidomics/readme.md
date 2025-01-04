# Class II elution experiments with Monash

## In this experiment, we pulled down I-Ab peptides from naive mice (3 mice pooled into 1)
| File                                | Vendor | Reference | Notes |
|-------------------------------------|----------|----------|----------|
| 202408_Pooled_Taconic_Bac_Query.csv | Taconic         |  See library below     |          |
| 202408_Pooled_Taconic_B6_Query.csv  | Taconic         |  _Mus musculus_           |          |
| 202408_Pooled_Jax_Bac_Query.csv     | Jackson         |  See library below         |          |
| 202408_Pooled_Jax_B6_Query.csv      | Jackson         | _Mus musculus_            |          |
| 202408_Pooled_Hill_Bac_Query.csv    | Hill         |  See library below         |          |
| 202408_Pooled_Hill_B6_Query.csv     | Hill         |  _Mus musculus_           |          |


## Results after running with NetMHCIIPan 4.3 (15 AA)
# Peptide Data Table
| Peptide                                 | Peptide Core  | WB/SB | %Rank | Taconic | Jackson | Hill | Taxa            | Protein / GenBank                                | Location |
|-----------------------------------------|---------------|-------|-------|---------|---------|------|-----------------|-------------------------------------------------|----------|
| MNYAWPTAEIAVM(+15.99)GG (oxidated Met)  | YAWPTAEIA     | SB    | 0.18  | Y       | N       | Y    | M. intestinale  | acyl-CoA carboxylase subunit beta / MEE0979697.1 | 421-435  |
| GDVSAYIPTNVISIT                         | YIPTNVISI     | WB    | 3.15  | Y       | Y       | Y    | D. dubosii / sp. B8                            |                                             |          |
| DVSAYIPTNVISITD                         | YIPTNVISI     | WB    | 1.37  | Y       | Y       | Y    | D. dubosii / sp. B8                            |                                             |          |
| VSAYIPTNVISITDG                         | YIPTNVISI     | WB    | 1.00  | Y       | Y       | Y    | D. dubosii / sp. B8                            | F0F1 ATP synthase subunit alpha / WP_136414477.1 | 349-364  |
| SAYIPTNVISITDGQ                         | YIPTNVISI     | WB    | 2.26  | Y       | N       | N    | P. intestinale                                  |                                             |          |
| LTLEAAPSVS(+79.97)EGVAR (phosphorylated Ser) | EAAPSVSEG   | SB    | 0.17  | N       | N       | Y    | P. intestinale                                  | UvrD-helicase domain-containing protein / WP_290374704.1 | 146-161  |
| TLEAAPSVSEGVARL                         | EAAPSVSEG     | SB    | 0.46  | N       | N       | Y    | P. intestinale                                  |                                             |          |
| LEAAPSVSEGVARLF                         | EAAPSVSEG     | WB    | 3.52  | N       | N       | Y    | P. intestinale                                  |                                             |          |
| KVRRIVNEPTAAALA                         | VNEPTAAAL     | WB    | 1.35  | N       | N       | Y    | M. intestinale                                  |                                             |          |
| VRRIVNEPTAAALAY                         | VNEPTAAAL     | SB    | 0.94  | N       | N       | Y    | M. intestinale                                  |                                             |          |
| RRIVNEPTAAALAYG                         | VNEPTAAAL     | SB    | 0.62  | N       | N       | Y    | M. intestinale                                  | molecular chaperone DnaK / WP_289823165.1     | 164-179  |
| RIVNEPTAAALAYGL                         | VNEPTAAAL     | WB    | 2.03  | N       | N       | Y    | M. intestinale                                  |                                             |          |



## Ideas on post-translationally modified peptides
One question is whether during oxidative stress, we see more PTM changes (e.g. oxidation, phosphorylation) that would general novel T-cell epitopes following allo SCT (vs. syn SCT).  It is widely known that a variety of autoimmune responses in human diseases depend on PTMs of self-proteins (e.g. citrullination in RA).  As another example with respect to T cell biology, the acetylation of the N-terminal peptide of myelin basic protein (MBP Ac1-11) is required for the development of experimental autoimmune encephalomyelitis (EAE), the murine model of human multiple sclerosis. 
 Here are some sources to look at:

| PMID                                | Year | Journal | Notes |
|-------------------------------------|----------|----------|----------|
| 36555449 | 2022         |  Int J Mol Sci.     | Emphasizes Citrullination, Homocitrullination and Acetylation, uses RA as disease model        |
| 11094420 | 2000         |  Arthritis Res.     | Focuses on autoantibodies       |
| 2431317 | 1986         |  Nature     | Encephalitogenic T cell clones specific for MBP Ac1-11 induce EAE, while the non-acetylated peptide fails to stimulate T cells or induce EAE      |




