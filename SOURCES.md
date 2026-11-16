# Published parameters and external data

The table distinguishes published parameter values implemented in the code from external data referenced in the study. Dataset links point to their original providers.

## Parameters used by the age model

| Model input | Public source and location |
| --- | --- |
| Atlantic bigeye growth and length dispersion | Waterhouse et al. (2022), *Updated growth models for bigeye tuna (Thunnus obesus) in the Atlantic Ocean*, Table 4 and Appendix 5 Table 1. The model uses the otolith-based von Bertalanffy growth parameters (161.21 cm, 0.392 yr⁻¹, −0.239 yr). The common 8.75 cm length SD approximates the 8.752 cm residual standard error from a separate fit with asymptotic length fixed to 160 cm and PFL age-1 and age-2 fish excluded. [NOAA full text](https://repository.library.noaa.gov/view/noaa/45309/noaa_45309_DS1.pdf) |
| Fork length to round weight | Parks et al. (1982), *Length–weight relations for bigeye tuna captured in the eastern Atlantic Ocean*, Table 7, relation 2, printed p. 223. The equation uses fork length in cm and round weight in kg. [IRD archive](https://horizon.documentation.ird.fr/exl-doc/pleins_textes/divers14-02/010028719.pdf) · [ICCAT adopted relationships, BET row](https://iccat.int/Documents/SCRS/Manual/Appendices/Appendix%204%20III%20Trop.pdf) |
| GG to RW | ICCAT (2021), *PSM Inspector Training Manual*, Appendix 1, Table 4, BET row. The factor for GG is 1.09. [ICCAT PDF](https://iccat.int/com2021/ENG/PWG_406_APP1_ENG.pdf) |
| Natural mortality | Lorenzen (2000), *Allometry of natural mortality as a basis for assessing optimal release size in fish-stocking programmes*, provides the inverse-length mortality form. The ICCAT 2025 Atlantic bigeye tuna data preparation report discusses average adult mortality of 0.32 yr⁻¹. The model calibrates the mean on monthly ages 4 through 10 yr. [Lorenzen paper](https://fisheriessolutions.org/wp-content/uploads/2015/02/Lorenzen_2000_Allometry-stocking.pdf) · [ICCAT report](https://www.iccat.int/Documents/Meetings/Docs/2025/Reports/2025_BET_DP_ENG.pdf) |
| External length-selectivity scenarios | Urtizberea et al. (2025), *Preliminary 2025 Stock Synthesis Assessment Model for Bigeye Tuna in the Atlantic Ocean*, Section 2.6 and Tables 2 and 5. Fleet 14 (tropical CTP longline) and fleet 15 (southern CTP longline) supply the two logistic parameter pairs. The logistic form is in Methot et al. (2020), *Stock Synthesis User Manual Version 3.30.15*, Eq. (21). [ICCAT assessment](https://www.iccat.int/Documents/CVSP/CV082_2025/n_5/CV082050151.pdf) · [NOAA manual](https://repository.library.noaa.gov/view/noaa/23884/noaa_23884_DS1.pdf) |

## External datasets cited for the companion study

| Dataset | Role in the study | Original source |
| --- | --- | --- |
| ICCAT Task 2 Size Samples (T2SZ), BET | Public sampled size frequencies used for an external size comparison. | [ICCAT Task 2 size](https://www.iccat.int/en/accesingdb.html) · [Field definitions](https://www.iccat.int/Data/t2size/ReadmeT2SZ-CS.pdf) |
| Andrews et al. (2024), *Age validation of yellowfin and bigeye tuna using post-peak bomb radiocarbon dating confirms long lifespans in the western and central Pacific Ocean* | Supplementary Table S4 provides Pacific bigeye age and length data used for an external comparison. | [DOI and supplementary material](https://doi.org/10.1093/icesjms/fsae074) |
| Global Ocean Ensemble Physics Reanalysis | Downstream environmental application. Product `GLOBAL_MULTIYEAR_PHY_ENS_001_031`, DOI `10.48670/moi-00024`. | [Copernicus Marine product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_ENS_001_031) |
| Global Ocean OSTIA Sea Surface Temperature and Sea Ice Reprocessed | Downstream environmental application. Product `SST_GLO_SST_L4_REP_OBSERVATIONS_010_011`, DOI `10.48670/moi-00168`. | [Copernicus Marine product page](https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011) |

The dataset citations are for provenance. The examples and the age-inversion package do not distribute or use these datasets.
