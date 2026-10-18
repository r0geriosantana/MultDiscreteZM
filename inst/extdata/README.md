# Data files

Data used in the applications of the paper (Section 5):

| File | Content |
|---|---|
| `Vertical_2015_Alagoas.xlsx` | 102 municipalities of Alagoas, 2015: `Municipio`, `Gestantes` (Y1, HIV-positive pregnant women), `Criancas` (Y2, children exposed to HIV), `Congenita` (Y3, congenital syphilis), `IDH` (municipal HDI), `GINI`. Maceió is excluded in the analyses (101 municipalities). |
| `dados_prontosMG.xlsx` | Minas Gerais, 76 municipalities, confirmed cases in 2023 (year of notification; year of diagnosis for leprosy) by municipality of residence, extracted from SINAN/DATASUS in 2025 (later extractions may differ slightly because SINAN records are revised): `Esquistossomose` (Y1, schistosomiasis), `LeishmanioseVisceral` (Y2, human visceral leishmaniasis), `Hanseniase` (Y3, leprosy), plus `cod_ibge`, `Populacao`, `GINI` and `IDH`. The paper scripts use the joint frequencies of the paper's table `freq_triv2`, embedded in `inst/paper/00_setup.R`. |
| `coorte_MG_2020.csv` | Classification of the 76 municipalities of the Minas Gerais cohort by their notifications in 2020 (SINAN/DATASUS, year of notification; year of diagnosis for leprosy; extracted in 2026): `cod_ibge`, `Municipio`, `Esquistossomose_2020`, `LeishmanioseVisceral_2020`, `Hanseniase_2020`, and `Grupo` = `Ay` (at least one notification in 2020, 57 municipalities) or `A0` (none, 19 municipalities). Used by `inst/paper/06_minas_selection.R`. |

Sources: Statistical Yearbook of the State of Alagoas 2015 (SESAU); DATASUS, Ministry of Health.
