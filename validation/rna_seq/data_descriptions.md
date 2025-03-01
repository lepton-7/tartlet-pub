# Transcriptionally-active Riboswitch Tracer Leveraging Edge deTection (TaRTLEt) validation datasets

### RNA-Seq data description

**e_coli_sastry2019_acc.txt**: Run accession list for *Escherichia coli* RNA-seq datasets from [Sastry 2019](https://www.nature.com/articles/s41467-019-13483-w). This list contains original runs from the paper (GSE122211, GSE122295, GSE122296, **but did not include GSE122320**) as well as datasets from [Guzman2018](https://www.biorxiv.org/content/10.1101/310946v2.full) (GSE114358) and [Anand2019](https://www.nature.com/articles/s41564-018-0340-2) (GSE122779).

**b_sub_guo2023_acc.txt**: Run accession list for *Bacillus subtilis* str. 168 from [Guo 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9879117/) ([GSE219221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219221)).

**b_sub_tilburg2022_acc.txt**: Run accession list for *Bacillus subtilis* str. 168 from [Tilburg 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9049611/) ([GSE169409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169409)). *These reads may not be paired since the fasterq-dump did not produce paired files*.

**b_sub_smith2023_acc.txt**: Run accession for *Bacillus subtilis* str. 168 from [GSE226559](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226559). List captures all runs 65.5GB.

---
*<u>Turns out HISAT-2 could not process these because they seem to be single-end reads packed as paired-end. Need to be thrown out</u>*

~~**p_fluor_wiesmann2022_acc.txt**: Run accession list for *P. fluorescens* from [Wiesmann 2022](https://www.nature.com/articles/s41396-022-01343-3) ([GSE190448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190448)).~~

~~**p_aeru_wiesmann2022_acc.txt**: Run accession list for *P. aeruginosa* from [Wiesmann 2022](https://www.nature.com/articles/s41396-022-01343-3) ([GSE190448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190448)).~~

---
**m_tuber_H37Rv_martini2023_acc.txt**: Run accession list for *M. tuberculosis* str H37Rv from [Martini 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10178195/) ([GSE218354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218354)).

**s_elon_lu2023_acc.txt**: Run accession list for *Synechococcus elongatus* from [Lu 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10365998/) ([GSE227397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227397)).

**c_vibrioides_mclaughlin2023_acc**: Run accession list for *Caulobacter vibrioides* from [McLaughlin 2023](https://pubmed.ncbi.nlm.nih.gov/37645952/) ([GSE241057](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241057)).

~~**a_fischeri_ES114_thompson2016_acc**: Run accession list for *Aliivibrio fischeri ES114* from [Thompson 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5409853/) ([GSE80607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80607)).~~ *These ended up being single ended reads.*

**a_fischeri_ES114_griend2023_acc**: Run accession list for *Aliivibrio fischeri ES114* from [Griend 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10441365/) ([GSE237189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237189)).

**b_theta_bolam2019_acc**: Run accession list for *Bacteroides thetaiotaomicron* from [Bolam 2019](https://www.nature.com/articles/s41564-019-0466-x) ([GSE129572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129572)).

**b_theta_lewis2023_acc**: Run accession list for *Bacteroides thetaiotaomicron* from [GSE169260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169260).

**s_spcc6803_wilde2022_wt_acc**: Run accession list for *Synechocystis sp. PCC 6803* from [GSE132275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132275). This list only captures wild-type runs.

**p_fluor_hasnain2023_acc**: Run accession list for *Pseudomonas fluorescens* from [Hasnain 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10229592/) ([GSE200822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200822)).

**s_coelicolor_wezel2023_wt_acc**: Run accesion list for *Streptomyces coelicolor* from [GSE234439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234439). This list only captures wild-type runs.

~~**a_rabiei_snrc2015_acc**: Run accession list for the fungus *Ascochyta rabiei* from [PRJNA288273](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA288273); [SRX1080581](https://www.ncbi.nlm.nih.gov/sra/SRX1080581[accn]). Of the 31 runs, only the ones in paired read layout are captured in the list. *Seems like the paired reads are split into individual SRA runs. The second in pair file was renamed according to first in pair's SRA ID*~~. *Did not work; not worth troubleshooting*.

**g_hirsutum_sun2023_acc**: Run accession list for C O T T O N *Gossypium hirsutum* from [GSE182982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182982).

**g_hirsutum_mei2022_acc**: Run accession list for C O T T O N *Gossypium hirsutum* from [GSE206663](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206663). >36Gbyte total; do not use unless the one above is insufficient.

**x_albi_wimmer2023_acc**: Run accession list for *Xanthomonas albilineans* from [GSE229478](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229478) 4.4GB.

**p_salmo_gabriela2023_acc**: Run accession list for *Piscirickettsia salmonis* from  [Carril2023](https://pubmed.ncbi.nlm.nih.gov/37953796/) ([GSE235725](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235725)). List captures all 6 single culture runs out of the total 9 runs 17.8GB.

**s_sanguinis_puccio2022_acc**: Run accession list for *Streptococcus sanguinis SK36* from [Puccio 2022](<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8844241/>) ([GSE174672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174672) 3.8GB).

**p_cholor_aureo3084_wang2016_acc**: Run accession list for *Pseudomonas chlororaphis* subsp. *aureofaciens* 30-84 from [Wang 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4727817/) ([GSE61200](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61200) 8.6GB).

**b_anth_corsi2020_acc**: Run accession list for *Bacillus anthracis* from [Corsi 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7843513/) ([GSE152356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152356)) 21.4GB.

~~**n_gonorr_zachary2021_acc**: Run accession list for *Neisseria gonorrhoeae* from [Zachary 2021](https://pubmed.ncbi.nlm.nih.gov/34515630/) ([GSE177032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177032) 5GB).~~ *These were single ended.*

**n_gonorr_cris2022_acc**: Run accession list for *Neisseria gonorrhoeae* from [Ray2022](https://journals.asm.org/doi/full/10.1128/iai.00009-22) ([GSE191020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191020)) 23.9GB.

**a_kunk_seeger2023_acc**: Run accession list for *Apilactobacillus kunkeei* from [Seeger 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10496945/) ([GSE205998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205998) 500MB).

**k_pneum_liu2023_acc**: Run accession list for *Klebsiella pneumoniae* from [GSE229867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229867) 2.1GB.

**s_enter_typh_kant2022_acc**: Run accession list for *Salmonella enterica* subsp. enterica serovar Typhimurium from [Kant2023](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002051) ([GSE203342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203342)) 13.9GB.

**a_baum_pokhrel2023_acc**: Run accession list for *Acinetobacter baumannii* from [GSE183334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183334) 3.8GB.

~~**c_basil_kugler2023_acc**: Run accession list for *Cupriavidus basilensis* from [GSE216827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216827). List without pgypsum and apatite condition 26.4GB.~~ *These were single ended.*

**m_smeg_grigorov2023_acc**: Run accession list for *Mycolicibacterium smegmatis* MC2 155 from [Grigorov 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10454040/) ([GSE232901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232901)). List captures only total RNA fraction 7.3GB.

**s_meli_fagorzi2021_acc**: Run accession list for *Sinorhizobium meliloti* from [Fagorzi 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7901481/) ([GSE151705](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151705)). List captures only the 11 no treatment conditions from the total 58 runs 9.7GB.

**d_vulg_gao2016_acc**: Run accession list for *Desulfovibrio vulgaris* str. Hildenborough from [Gao 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5007762/) ([GSE78834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78834)) 2.8GB.

**d_vulg_guo2019_acc**: Run accession list for *Desulfovibrio vulgaris* str. Hildenborough from [Chen2019](https://www.sciencedirect.com/science/article/pii/S0160412018326138) ([GSE101911](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101911)) 7.9GB.

**b_frag_fiebig2023_acc**: Run accession list for *Bacteroides fragilis* from [Fiebig 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10197588/) ([GSE220692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220692)) 3.1GB.

**e_fae_avican2021_acc**: Run accession list for *Enterococcus faecalis* from [Avican 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8172932/) ([GSE152295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152295)). List captures 2/3 replicates for each of the non-untreated conditions and all of the untreated conditions 15.9GB.

**b_xyla_despres2016_acc**: Run accession list for *Bacteroides xylanisolvens* from [Despres 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4769552/) ([GSE74379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74379)) 23.9Gb.

**x_ory_park2022_acc**: Run accession list for *Xanthomonas oryzae* pv. oryzae KACC 10331 from [Park 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9709210/) ([GSE89651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89651)) 11.4GB.

**s_epi_avican2021_acc**: Run accession list for *Staphylococcus epidermidis* from [Avican 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8172932/) ([GSE152295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152295)) 23.5GB.

**b_pseudo_avican2021_acc**: Run accession list for *Burkholderia pseudomallei* from [Avican 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8172932/) ([GSE152295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152295)) 23.5GB.

**c_diff_weiss2021_acc**: Run accession list for *Clostridioides difficile* from [Weiss 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8767335/) ([GSE173804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173804)). List only captures the 6 runs of strain CD630 19.8GB.

**e_limo_jeong2020_acc**: Run accession list for *Eubacterium limosum* from [GSE149269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149269) 6.2GB.

## Validation Survey of *B. subtilis* and *E. coli*

### Alignment

Indexes were generated via:

```bash
hisat2-build -p 48 validation/genomes/Bsubtilis168.GCF_000009045.1.fna validation/survey/b_sub_168_index/b_sub_168_index
```

## Miscellaneous

Start an interactive shell on pitzer:

```bash
srun -n 1 --cpus-per-task=48 --account=PDS0325 --time=1:00:00 --pty /bin/bash
```

For updated list of validation set:

```python
import pandas as pd
from pathlib import Path

df = pd.read_csv(Path("validation/tables/inf_results.csv"))
for x in sorted(pd.unique(df["dataset"])):
    print(x)
```

Move job files to their own directories:

```bash
DSETS=(
    # "a_baum"
    "a_fischeri_ES114"
    "a_kunk"
    "b_anth"
    "b_frag"
    "b_pseudo"
    "b_sub_168"
    "b_theta"
    "b_xyla"
    # "c_basil"
    "c_diff"
    "c_vibrioides"
    "d_vulg"
    "e_coli"
    "e_fae"
    "e_limo"
    "k_pneum"
    "m_smeg"
    "m_tuber"
    "n_gonorr"
    "p_aeru"
    "p_cholor_aureo3084"
    "p_fluor"
    "p_salmo"
    "s_coelicolor"
    "s_elon"
    "s_enter_typh"
    "s_epi"
    "s_meli"
    "s_sanguinis"
    "s_spcc6803"
    "x_albi"
    "x_ory"
)

for DSET in ${DSETS[@]}; do
    mkdir -p "./jobs/align/$DSET"
    mv ./jobs/*$DSET.* ./jobs/align/$DSET
done
```
