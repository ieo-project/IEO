---
output:
  BiocStyle::html_document
---

<!---
The following chunk of code, which should not be shown in the resulting document (echo=FALSE)
sets up global processing options, such as forcing 'knitr' to stop when an error
in the R code is encountered, caching of the results in the 'cache'
directory and asking 'knitr' to figure out automatically the dependencies among
code chunks to re-calculate cached results (autodep=TRUE).

Other options could be changing the name of the directory where figures end up
('figure' by default), etc. For a full account of 'knitr' options please consult
http://yihui.name/knitr/options

At the end of the chunk a 'cat()' call is made to dump a CSS file that gives
a better look-and-feel than the knitr default one. See the source css/ieo.css
and the resulting projectTemplate.html to understand where this is being dumpted.
--->




# Quality assessment

## Data import

We start importing the raw table of counts.

<!--
The option 'message=FALSE' avoid dumping R messages such as "Loading required package: methods"
into the output of the report.
-->



```r
library(SummarizedExperiment)
setwd("/home/guille/Guille/MBHS/3T/IEO/project")
se <- readRDS("seKIRC.rds")
# se
```

```
class: RangedSummarizedExperiment 
dim: 20115 614 
metadata(5): experimentData annotation cancerTypeCode cancerTypeDescription
  objectCreationDate
assays(1): counts
rownames(20115): 1 2 ... 102724473 103091865
rowData names(3): symbol txlen txgc
colnames(614): TCGA.3Z.A93Z.01A.11R.A37O.07 TCGA.6D.AA2E.01A.11R.A37O.07 ...
  TCGA.CZ.5988.11A.01R.1672.07 TCGA.CZ.5989.11A.01R.1672.07
colData names(549): type bcr_patient_uuid ... lymph_nodes_aortic_pos_by_ihc
  lymph_nodes_aortic_pos_total
```

Explore the column (phenotypic) data, which in this case corresponds to clinical
variables, and their corresponding metadata.


```r
dim(colData(se))
```

```
[1] 614 549
```

```r
colData(se)[1:5, 1:5]
```

```
DataFrame with 5 rows and 5 columns
                                 type                     bcr_patient_uuid bcr_patient_barcode
                             <factor>                             <factor>            <factor>
TCGA.3Z.A93Z.01A.11R.A37O.07    tumor 2B1DEA0A-6D55-4FDD-9C1C-0D9FBE03BD78        TCGA-3Z-A93Z
TCGA.6D.AA2E.01A.11R.A37O.07    tumor D3B47E53-6F40-4FC8-B5A4-CBE548A770A9        TCGA-6D-AA2E
TCGA.A3.3306.01A.01R.0864.07    tumor 9fb55e0b-43d8-40a3-8ef2-d198e6290551        TCGA-A3-3306
TCGA.A3.3307.01A.01R.0864.07    tumor 7ac1d6c6-9ade-49af-8794-10b5b96b2b05        TCGA-A3-3307
TCGA.A3.3308.01A.02R.1325.07    tumor 3cbca837-f5a7-4a87-8f02-c59eac232d5a        TCGA-A3-3308
                             form_completion_date prospective_collection
                                         <factor>               <factor>
TCGA.3Z.A93Z.01A.11R.A37O.07           2014-11-11                    YES
TCGA.6D.AA2E.01A.11R.A37O.07            2014-3-17                    YES
TCGA.A3.3306.01A.01R.0864.07            2010-8-23                     NO
TCGA.A3.3307.01A.01R.0864.07            2010-4-13                     NO
TCGA.A3.3308.01A.02R.1325.07            2010-4-12                     NO

```

```r
mcols(colData(se), use.names=TRUE)
colnames(colData(se))

## Filter by NA ##
df <- data.frame(colData(se))

NA_count <- colSums(df == "[Not Available]" | is.na(df) | df == "[Not Applicable]")
NA_count[(NA_count <300)]

colSums(is.na(df))
na_count <- colSums(is.na(df))
na_count[(na_count <300)]

# example to obtain CDEID
#which(colnames(colData(se)) == "retrospective_collection")

summary(colData(se)[,which(colnames(colData(se)) == "initial_pathologic_dx_year")])
mcols(colData(se), use.names=TRUE)[which(colnames(colData(se)) == "initial_pathologic_dx_year"),]


```

```
DataFrame with 549 rows and 2 columns
                                                         labelDescription       CDEID
                                                              <character> <character>
type                                           sample type (tumor/normal)          NA
bcr_patient_uuid                                         bcr patient uuid          NA
bcr_patient_barcode                                   bcr patient barcode     2673794
form_completion_date                                 form completion date          NA
prospective_collection            tissue prospective collection indicator     3088492
...                                                                   ...         ...
lymph_nodes_pelvic_pos_total                               total pelv lnp     3151828
lymph_nodes_aortic_examined_count                           total aor lnr     3104460
lymph_nodes_aortic_pos_by_he                          aln pos light micro     3151832
lymph_nodes_aortic_pos_by_ihc                                 aln pos ihc     3151831
lymph_nodes_aortic_pos_total                                total aor-lnp     3151827
```

These metadata consists of two columns of information about the clinical variables.
One called `labelDescription` contains a succint description of the variable, often
not more self-explanatory than the variable name itself, and the other called
'CDEID' corresponds to the so-called `Common Data Element (CDE)` identifier. This
identifier can be use in https://cdebrowser.nci.nih.gov to search for further
information about the associated clinical variable using the `Advanced search`
form and the `Public ID` attribute search.

Now, explore the row (feature) data.


```r
rowData(se)
```

```
DataFrame with 20115 rows and 3 columns
               symbol     txlen              txgc
          <character> <integer>         <numeric>
1                A1BG      3322 0.564419024683925
2                 A2M      4844 0.488232865400495
9                NAT1      2280 0.394298245614035
10               NAT2      1322 0.389561270801815
12           SERPINA3      3067 0.524942940984676
...               ...       ...               ...
100996331       POTEB      1706 0.430832356389215
101340251    SNORD124       104 0.490384615384615
101340252   SNORD121B        81 0.407407407407407
102724473      GAGE10       538 0.505576208178439
103091865   BRWD1-IT2      1028 0.592412451361868
```

```r
rowRanges(se)
```

```
GRanges object with 20115 ranges and 3 metadata columns:
            seqnames            ranges strand |      symbol     txlen              txgc
               <Rle>         <IRanges>  <Rle> | <character> <integer>         <numeric>
          1    chr19 58345178-58362751      - |        A1BG      3322 0.564419024683925
          2    chr12   9067664-9116229      - |         A2M      4844 0.488232865400495
          9     chr8 18170477-18223689      + |        NAT1      2280 0.394298245614035
         10     chr8 18391245-18401218      + |        NAT2      1322 0.389561270801815
         12    chr14 94592058-94624646      + |    SERPINA3      3067 0.524942940984676
        ...      ...               ...    ... .         ...       ...               ...
  100996331    chr15 20835372-21877298      - |       POTEB      1706 0.430832356389215
  101340251    chr17 40027542-40027645      - |    SNORD124       104 0.490384615384615
  101340252     chr9 33934296-33934376      - |   SNORD121B        81 0.407407407407407
  102724473     chrX 49303669-49319844      + |      GAGE10       538 0.505576208178439
  103091865    chr21 39313935-39314962      + |   BRWD1-IT2      1028 0.592412451361868
  -------
  seqinfo: 455 sequences (1 circular) from hg38 genome
```

To perform quality assessment and normalization we need first to load the
[edgeR](http://bioconductor.org/packages/edgeR) R/Bioconductor package and
create a `DGEList` object.


```r
library(edgeR)

dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
head(dge$samples)
```

```
Warning in as.data.frame(x, row.names = NULL, optional = optional, ...):
Arguments in '...' ignored
```

```r
saveRDS(dge, "dge.rds")
```

Now calculate $\log_2$ CPM values of expression and put them as an additional
assay element to ease their manipulation.


```r
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM[1:5, 1:5]
```

```
   TCGA.3Z.A93Z.01A.11R.A37O.07 TCGA.6D.AA2E.01A.11R.A37O.07 TCGA.A3.3306.01A.01R.0864.07
1                      3.681094                     1.308506                   -0.2268151
2                     11.070451                     9.467119                    9.8676287
9                     -6.544714                    -6.544714                   -6.5447145
10                    -6.544714                    -6.544714                   -6.5447145
12                     5.196098                     5.386859                    4.6330031
   TCGA.A3.3307.01A.01R.0864.07 TCGA.A3.3308.01A.02R.1325.07
1                      0.203261                    0.2380578
2                     11.166993                   11.2922910
9                     -6.544714                   -6.5447145
10                    -6.544714                   -6.5447145
12                     5.479129                    5.8638782
```

## Sequencing depth

Let's examine the sequencing depth in terms of total number of sequence read counts
mapped to the genome per sample. Figure \@ref(fig:libsizes) below shows the
sequencing depth per sample, also known as library sizes, in increasing order.

<!---
you can control the height and width in pixels of the figure with 'out.height' and
'out.width'. Figures are automatically numbered, to refer to them in the main test
you should use the notation shown above as \@ref(fig:xxxx) with xxxx being the label
in the code chunk that also gives the filename of the figure. This name must be unique
--->

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/libsizes-1.png" alt="Library sizes in increasing order." width="600px" />
<p class="caption">(\#fig:libsizes)Library sizes in increasing order.</p>
</div>
This figure reveals substantial differences in sequencing depth between samples
and we may consider discarding those samples whose depth is substantially lower
than the rest. To identify who are these samples we may simply look at the
actual numbers including portion of the sample identifier that distinguishes them.


```r
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
names(sampledepth) <- substr(colnames(se), 6, 12)
sort(sampledepth)

## what to do with NA's?????? ##

ord <- order(dge$sample$lib.size)
dge$sample$lib.size
barplot(dge$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples", col = c("red", "blue")[se$gender[ord]], border=NA, space=c(0,0))
legend("topleft", c("female", "male", "NA"), fill = c("red", "blue"), inset = 0.02)
```

```
CJ.4642 AK.3454 B0.4697 B2.3923 B0.4841 B0.5706 B8.A54D A3.3326 B0.5702 B2.3924 B2.5633 
    3.0    14.4    14.5    15.1    15.6    17.2    17.7    18.1    18.3    18.9    19.4 
B2.3924 AK.3453 G6.A8L7 B8.A54F B2.5635 BP.5185 BP.4334 AK.3426 BP.4335 B0.4707 AK.3428 
   19.7    20.0    20.1    20.3    20.7    20.8    21.0    24.1    24.4    25.1    25.2 
CZ.5987 B2.5639 B0.4833 BP.4776 B0.4718 CJ.4638 A3.3316 G6.A8L6 B8.4619 BP.4992 CZ.5985 
   25.8    25.9    26.3    26.3    26.4    26.6    26.9    27.2    27.2    27.4    27.5 
B0.4813 BP.5187 AK.3433 3Z.A93Z B2.4098 B8.4620 B0.4847 B0.4843 BP.4160 MM.A564 B0.5117 
   27.6    27.6    27.7    27.8    27.8    28.0    28.3    28.4    28.4    28.5    28.6 
AK.3460 B0.4817 CJ.4873 B0.4828 B8.A54I G6.A5PC BP.5010 CJ.4639 B0.4849 BP.5191 B0.4822 
   28.8    29.1    29.5    29.7    30.0    30.0    30.3    30.4    30.5    30.5    30.6 
BP.4164 B0.5694 B0.4713 B0.5701 DV.A4VX B4.5832 BP.4804 B0.4700 CJ.4882 AK.3434 CZ.5987 
   30.6    30.7    30.8    30.9    30.9    31.0    31.0    31.1    31.2    31.4    31.4 
CJ.5672 CZ.5986 EU.5905 MM.A84U AK.3447 CJ.4912 6D.AA2E CJ.4875 B0.5691 BP.4159 B0.4714 
   31.5    31.5    31.5    31.6    31.7    31.7    32.0    32.1    32.4    32.5    32.8 
MM.A563 CJ.4636 CZ.5988 B8.A7U6 CJ.4878 A3.3306 B8.A54E BP.4163 CZ.5989 CZ.5988 AK.3431 
   32.8    32.9    32.9    33.0    33.0    33.2    33.2    33.2    33.2    33.2    33.3 
B0.4824 AK.3451 CJ.4635 B8.4153 CZ.5456 B0.4842 EU.5907 A3.3322 B0.5083 B8.4622 AK.3456 
   33.3    33.5    33.5    33.6    33.6    33.8    33.9    34.0    34.1    34.2    34.3 
B0.4819 B0.5703 DV.5566 BP.4354 CJ.5681 AK.3458 B0.5085 B0.5697 AK.3429 B0.5700 BP.4165 
   34.5    34.5    34.5    34.6    34.6    34.7    34.7    34.8    35.0    35.1    35.1 
BP.4167 CJ.5676 CZ.5986 B0.4811 B0.4844 BP.5169 A3.3323 BP.4799 B0.5092 A3.3383 B2.3923 
   35.1    35.1    35.2    35.3    35.4    35.4    35.5    35.5    35.6    35.7    35.7 
CJ.5682 AK.3461 A3.3307 CJ.4900 B8.5546 B8.A54G B0.4693 A3.3362 B0.4845 A3.3358 CZ.4862 
   35.7    35.8    35.9    35.9    36.0    36.2    36.4    36.5    36.6    36.6    36.7 
G6.A8L8 CJ.5689 B8.4143 DV.A4W0 CJ.5679 B0.4945 B8.4619 B8.5550 B0.4848 BP.5001 B0.4696 
   36.7    36.7    36.8    36.8    36.9    37.1    37.1    37.1    37.3    37.3    37.4 
B0.4834 B2.3923 CJ.6033 A3.3359 B2.5633 BP.4981 B0.5699 CW.5581 CJ.5677 DV.A4VZ BP.4807 
   37.4    37.4    37.4    37.5    37.5    37.5    37.6    37.6    37.7    37.9    38.1 
CJ.6028 CZ.5466 BP.4158 B0.5706 CZ.5459 A3.3308 BP.4162 BP.5184 AK.3445 B8.4151 CZ.4857 
   38.1    38.1    38.2    38.2    38.3    38.5    38.5    38.5    38.6    38.7    38.7 
AK.3425 CJ.5678 B8.A8YJ CJ.4868 A3.A8OV BP.5168 B2.5641 BP.4166 BP.4790 B0.5703 A3.3351 
   38.8    38.8    38.9    38.9    39.0    39.0    39.0    39.1    39.1    39.2    39.3 
B0.4836 B0.5707 A3.A8CQ B8.4622 BP.4988 AK.3455 CZ.5451 B0.4839 B0.5096 A3.3374 A3.A6NL 
   39.3    39.3    39.4    39.8    39.8    40.0    40.0    40.1    40.1    40.2    40.2 
BP.5170 A3.3387 B0.4701 B8.A54K BP.4327 B4.5843 BP.5190 DV.5565 CJ.4637 A3.3365 BP.5175 
   40.4    40.4    40.5    40.5    40.5    40.6    40.7    40.7    40.8    40.9    41.0 
BP.4337 BP.4351 A3.3347 AK.3427 B0.4814 B0.4838 B2.4101 AK.3465 CJ.4887 B0.4837 B0.5106 
   41.1    41.1    41.2    41.2    41.3    41.3    41.3    41.4    41.5    41.6    41.6 
MW.A4EC A3.3349 BP.4341 B0.5109 BP.4343 CJ.5683 BP.4346 CJ.4895 CJ.5678 B4.5838 BP.5177 
   41.8    41.9    41.9    42.0    42.0    42.0    42.1    42.1    42.1    42.2    42.2 
B0.5107 B8.4146 CJ.4870 CJ.4901 A3.3313 BP.5183 BP.4355 B2.A4SR BP.4960 CJ.4876 CJ.4920 
   42.3    42.3    42.3    42.3    42.4    42.4    42.5    42.6    42.6    42.7    42.7 
A3.A8OW B0.4691 B0.5121 B2.4099 CJ.4641 CJ.4881 BP.4759 DV.5573 B0.4700 B0.5692 A3.A8OU 
   42.8    42.9    42.9    42.9    42.9    42.9    43.0    43.0    43.1    43.1    43.2 
CZ.5982 A3.3329 B2.5633 B0.5099 B2.5635 B0.5693 BP.4344 AK.3436 BP.4325 BP.4787 CZ.4858 
   43.2    43.3    43.3    43.4    43.4    43.5    43.5    43.7    43.7    43.7    43.7 
CW.5584 A3.A6NJ B8.4148 BP.4993 DV.5575 B0.5696 B0.5110 B8.5551 CJ.4644 CW.6097 A3.A6NN 
   43.8    44.1    44.1    44.1    44.1    44.2    44.3    44.3    44.4    44.4    44.5 
B0.5095 BP.4798 B2.3924 B8.4620 CJ.4891 B0.5097 BP.4763 CJ.5689 BP.4332 BP.4969 A3.3317 
   44.5    44.6    44.8    44.9    45.1    45.2    45.2    45.2    45.3    45.3    45.4 
B0.5100 BP.4756 BP.4777 CJ.5684 B0.5690 BP.4784 BP.5000 CJ.4634 CW.5588 B0.4706 BP.4174 
   45.4    45.4    45.4    45.5    45.6    45.6    45.6    45.6    45.6    45.7    45.7 
AK.3440 CJ.5671 DV.5568 B0.4821 CJ.4923 BP.4769 A3.3382 B2.4102 CJ.4886 CW.5591 B0.5694 
   45.8    45.8    45.8    46.0    46.0    46.1    46.2    46.2    46.2    46.2    46.2 
BP.5174 CJ.4872 BP.4349 BP.4353 AK.3450 B4.5844 BP.4161 BP.4771 A3.3335 BP.4781 CJ.5680 
   46.3    46.3    46.5    46.5    46.6    46.6    46.6    46.6    46.7    46.7    46.7 
B0.5116 B2.5635 B8.A54J BP.5196 B0.4710 B2.5636 B8.5549 BP.4761 DV.A4W0 CZ.5984 A3.3320 
   46.8    46.8    46.8    46.8    46.9    46.9    46.9    46.9    47.0    47.0    47.1 
B0.5120 B0.5710 BP.4803 A3.3324 A3.3343 BP.4963 CZ.5469 CJ.5677 BP.4989 B0.5119 BP.4774 
   47.1    47.1    47.1    47.3    47.3    47.4    47.4    47.8    47.9    48.1    48.1 
CW.5589 B2.5636 BP.4345 CZ.4865 BP.4176 CJ.4874 B0.5080 BP.4971 CJ.4892 CZ.5452 A3.3378 
   48.2    48.2    48.3    48.3    48.5    48.5    48.6    48.6    48.6    48.6    48.7 
B0.4815 B0.5399 BP.4352 DV.5569 CZ.5989 BP.4758 B0.5115 BP.4782 BP.5178 CJ.4918 CJ.6027 
   48.7    48.7    48.7    48.7    48.7    48.8    48.9    48.9    48.9    48.9    48.9 
B0.4688 CJ.5680 B0.4846 B0.5400 CZ.5982 B0.5712 CJ.4643 CZ.4854 B0.5701 B0.5084 B0.5088 
   49.0    49.0    49.1    49.1    49.1    49.2    49.2    49.2    49.2    49.3    49.3 
BP.5181 B0.5711 CJ.4889 B0.4823 BP.4775 CJ.4884 B0.4712 BP.4801 CJ.4871 A3.A6NI B0.4703 
   49.3    49.5    49.5    49.6    49.6    49.6    49.7    49.7    49.7    49.8    49.8 
B8.5553 AS.3778 B0.4690 B0.5705 BP.4762 CW.6090 A3.3311 A3.3346 A3.3372 CJ.5681 BP.5186 
   49.8    49.9    49.9    49.9    50.0    50.0    50.1    50.1    50.1    50.1    50.2 
B8.5162 BP.4329 B0.5690 A3.3380 B8.4154 BP.5009 B0.5102 B8.5545 CW.5590 CZ.4861 B0.5094 
   50.3    50.3    50.3    50.5    50.5    50.5    50.6    50.6    50.6    50.6    50.7 
BP.4760 CJ.4885 CW.5584 B0.5699 CZ.5460 B0.4818 CJ.4894 CW.5589 B0.5108 DV.5576 EU.5904 
   50.7    50.7    50.7    50.8    50.8    50.9    50.9    50.9    51.0    51.0    51.0 
B0.5713 CJ.6031 AS.3777 CJ.4640 BP.4326 BP.5198 B2.5641 B8.5549 BP.5173 BP.5201 CW.6093 
   51.1    51.1    51.2    51.2    51.3    51.3    51.4    51.4    51.5    51.5    51.5 
CW.5591 B0.5113 CW.6087 CJ.4916 B0.5712 CZ.5468 CZ.5456 A3.3331 B0.4694 BP.4173 CJ.5679 
   51.5    51.6    51.6    51.7    51.9    52.0    52.1    52.2    52.2    52.2    52.2 
B0.5104 BP.4331 BP.5180 CW.5587 A3.A8OX BP.4789 CZ.5452 BP.4768 DV.5574 B0.5698 CJ.4888 
   52.4    52.4    52.4    52.4    52.5    52.7    52.7    52.8    52.9    53.0    53.0 
CZ.5453 A3.3385 CJ.5675 B0.4810 BP.4962 BP.4330 BP.4987 B8.5552 B4.5836 B0.5696 B0.5812 
   53.0    53.1    53.1    53.2    53.3    53.4    53.4    53.5    53.6    53.7    53.8 
B0.5709 B0.5081 CJ.6030 CZ.4860 BP.5195 A3.3363 CZ.5468 B0.5098 B8.4621 CJ.4890 GK.A6C7 
   53.9    54.0    54.0    54.0    54.1    54.3    54.4    54.6    54.6    54.7    54.7 
A3.3387 CW.5581 BP.4177 BP.4342 AK.3443 BP.5192 CZ.5467 BP.4765 BP.5006 B8.5552 BP.5007 
   54.8    54.8    54.9    55.0    55.2    55.2    55.3    55.4    55.5    55.6    55.6 
CW.5587 B0.5402 BP.4170 CJ.5672 CZ.5985 T7.A92I CW.5585 CZ.4863 CZ.4864 CZ.5467 A3.3319 
   55.6    55.6    55.7    55.8    56.1    56.1    56.1    56.2    56.2    56.2    56.3 
A3.3357 CJ.5676 CZ.5458 BP.4795 B0.5077 BP.5176 BP.5199 CJ.4903 A3.3358 BP.4766 DV.5567 
   56.3    56.3    56.4    56.5    56.6    56.6    56.6    56.6    56.7    56.7    56.7 
BP.4340 A3.3328 CJ.4908 CJ.6033 CZ.5464 CJ.4869 BP.4983 CZ.5984 BP.4964 CJ.4893 CW.6090 
   56.8    56.9    56.9    56.9    57.1    57.2    57.3    57.5    57.6    57.7    57.8 
B4.5378 B8.5165 B0.4816 CZ.5457 B0.5697 BP.4968 B0.5711 A3.3367 BP.4797 CJ.5686 BP.4967 
   57.9    57.9    58.0    58.1    58.2    58.2    58.2    58.3    58.3    58.4    58.6 
BP.4977 BP.4985 CW.5583 B0.4698 CZ.5470 CJ.4907 B8.A54H B0.5709 B0.5695 B4.5834 BP.5182 
   58.6    58.6    58.6    58.7    58.7    58.8    58.9    59.0    59.1    59.1    59.1 
CW.5580 B0.4827 BP.4959 B0.5075 CW.5585 B8.5158 B0.5705 CZ.5457 CJ.6032 BP.4986 CZ.4864 
   59.5    59.6    59.7    60.0    60.0    60.4    60.5    60.5    60.6    60.8    60.8 
CZ.5455 CZ.5462 CZ.5451 CJ.4905 B8.5164 BP.5008 CW.6088 A3.3373 BP.5202 CZ.4853 B0.5691 
   60.9    61.0    61.2    61.3    61.4    61.4    61.4    61.7    61.8    62.1    62.1 
CZ.5463 B0.5402 CJ.6030 B4.5835 BP.4973 B8.5159 CZ.5453 B0.4712 BP.4347 BP.4998 BP.4338 
   62.1    62.6    62.7    62.9    62.9    63.0    63.1    63.3    63.4    63.4    63.7 
CZ.4856 BP.5194 CZ.5465 CZ.5469 CJ.4902 BP.4770 BP.4982 CJ.4899 A3.3370 BP.4169 CZ.5461 
   63.8    64.0    64.1    64.3    64.4    64.5    64.6    64.6    64.7    65.0    65.4 
BP.4991 CW.5580 CZ.5470 BP.4961 A3.3325 BP.5004 CW.6087 BP.4999 CZ.5454 B4.5377 CZ.4859 
   65.7    65.8    66.0    66.3    66.5    66.6    66.6    66.8    67.2    67.4    67.5 
CJ.4904 CZ.4863 BP.4995 CW.6088 BP.4994 BP.5200 CZ.5466 A3.3352 CZ.5461 BP.4974 CJ.4897 
   67.6    67.6    67.7    67.8    68.2    69.0    69.0    69.3    69.5    69.7    72.1 
EU.5906 BP.4976 B0.4852 A3.3376 B8.5163 CZ.5465 BP.4975 BP.4970 CZ.5455 BP.4972 B0.4699 
   72.1    72.4    72.6    72.8    73.9    74.1    74.4    76.1    76.1    77.1    77.7 
CZ.4865 CZ.5463 BP.5189 BP.4965 CZ.4866 CZ.5462 CZ.5454 CZ.5458 AK.3444 
   77.7    78.2    78.4    78.7    79.4    80.1    80.6    81.9   115.5 
```

## Distribution of expression levels among samples

Let's look at the distribution of expression values per sample in terms of
logarithmic CPM units. Due to the large number of samples, we display tumor
and normal samples separately, and are shown in Figure \@ref(fig:distRawExp)

<!---
the option echo=FALSE hides the R code. When plotting in general one does not
want to see the code. Options fig.height and fig.width control height and width
of the plot in inches while out.height and out.width do it in the final output
file; see http://yihui.name/knitr/options for full details.
--->

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/distRawExp-1.png" alt="Non-parametric density distribution of expression profiles per sample." width="800px" />
<p class="caption">(\#fig:distRawExp)Non-parametric density distribution of expression profiles per sample.</p>
</div>

We do not appreciate substantial differences between the samples in the
distribution of expression values.

## Distribution of expression levels among genes

Let's calculate now the average expression per gene through all the samples.
Figure \@ref(fig:exprdist) shows the distribution of those values across genes.

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/exprdist-1.png" alt="Distribution of average expression level per gene." width="400px" />
<p class="caption">(\#fig:exprdist)Distribution of average expression level per gene.</p>
</div>

## Filtering of lowly-expressed genes

In the light of this plot, we may consider a cutoff of 1 log CPM unit as minimum value
of expression to select genes being expressed across samples. Using this cutoff we proceed
to filter out lowly-expressed genes.


```r
mask <- avgexp > 1
dim(se)
```

```
[1] 20115    91
```

```r
se.filt <- se[mask, ]
dim(se.filt)
```

```
[1] 11368    91
```

```r
dge.filt <- dge[mask, ]
dim(dge.filt)
```

```
[1] 11368    91
```

Store un-normalized versions of the filtered expression data.


```r
saveRDS(se.filt, file.path("results", "se.filt.unnorm.rds"))
saveRDS(dge.filt, file.path("results", "dge.filt.unnorm.rds"))
```

## Normalization

We calculate now the normalization factors on the filtered expression data set.


```r
dge.filt <- calcNormFactors(dge.filt)
```

Replace the raw log2 CPM units in the corresponding assay element of the `SummarizedExperiment`
object, by the normalized ones.


```r
assays(se.filt)$logCPM <- cpm(dge.filt, log=TRUE, normalized.lib.sizes=TRUE, prior.count=0.25)
```

Store normalized versions of the filtered expression data.


```r
saveRDS(se.filt, file.path("results", "se.filt.rds"))
saveRDS(dge.filt, file.path("results", "dge.filt.rds"))
```

## MA-plots

We examine now the MA-plots of the normalized expression profiles. We look first to
the tumor samples in Figure \@ref(fig:maPlotsTumor).

<!---
Here we make a MA-plot for each sample. The options 'fig.height' and 'fig.width'
control the relative image size in *inches*. The final image size results from
'height'x'dpi' and 'width'x'dpi', where 'dpi' is the image resolution in
"dots per inch" (by default dpi=72). To scale the image to a desired size use
'out.width' and 'out.height'. More information at http://yihui.name/knitr/options
--->

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/maPlotsTumor-1.png" alt="MA-plots of the tumor samples." width="600" />
<p class="caption">(\#fig:maPlotsTumor)MA-plots of the tumor samples.</p>
</div>

We do not observe samples with major expression-level dependent biases. Let's
look now to the normal samples in Figure \@ref(fig:maPlotsNormal).

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/maPlotsNormal-1.png" alt="MA-plots of the normal samples." width="600" />
<p class="caption">(\#fig:maPlotsNormal)MA-plots of the normal samples.</p>
</div>

We do not observe either important expression-level dependent biases among the normal samples.

## Batch identification

We will search now for potential surrogate of batch effect indicators. Given that each sample
names corresponds to a TCGA barcode (see https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode),
following the strategy described in http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview
we are going to derive different elements of the TCGA barcode and examine their distribution
across samples.


```r
tss <- substr(colnames(se.filt), 6, 7)
table(tss)
```

```
tss
KL KM KN KO 
30  9 36 16 
```

```r
center <- substr(colnames(se.filt), 27, 28)
table(center)
```

```
center
07 
91 
```

```r
plate <- substr(colnames(se.filt), 22, 25)
table(plate)
```

```
plate
2315 2403 
  90    1 
```

```r
portionanalyte <- substr(colnames(se.filt), 18, 20)
table(portionanalyte)
```

```
portionanalyte
01R 11R 21R 
 25  65   1 
```

```r
samplevial <- substr(colnames(se.filt), 14, 16)
table(samplevial)
```

```
samplevial
01A 11A 
 66  25 
```

From this information we can make the following observations:

  * All samples were sequenced at the same center

  * All samples belong to one of two combinations of tissue type and vial, matching the
    expected tumor and normal distribution.

  * Samples were collected across different tissue source sites (TSS).

  * All samples were sequenced within the same plate, except for the following one:


```r
colnames(se.filt)[plate == "2403"]
```

```
[1] "TCGA.KM.8639.01A.11R.2403.07"
```

  * All samples were sequenced using one of two portion and analyte combinations except fo the
    following one:


```r
colnames(se.filt)[portionanalyte == "21R"]
```

```
[1] "TCGA.KL.8323.01A.21R.2315.07"
```

We are going to use the TSS as surrogate of batch effect indicator. Considering our outcome
of interest as molecular changes between sample types, tumor vs. normal, we will examine now
the cross-classification of this outcome with TSS.


```r
table(data.frame(TYPE=se.filt$type, TSS=tss))
```

```
        TSS
TYPE     KL KM KN KO
  normal  6  0 17  2
  tumor  24  9 19 14
```

Observe that normal tissues with `TSS=KM` or `TSS=KO` are under-represented with respect to
the tumor tissues. If TSS is a source of expression variability, this under-representation
of those two TSS in the normal samples may lead to a potential confounding effect.

We examine now how samples group together by hierarchical clustering and multidimensional
scaling, annotating the outcome of interest and the the surrogate of batch indicator. We
calculate again log CPM values with a higher prior count to moderate extreme fold-changes
produced by low counts. The resulting dendrogram is shown in Figure \@ref(fig:sampleClustering).


```r
logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se.filt)
outcome <- paste(substr(colnames(se.filt), 9, 12), as.character(se.filt$type), sep="-")
names(outcome) <- colnames(se.filt)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))
```

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/sampleClustering-1.png" alt="Figure S6: Hierarchical clustering of the samples." width="1400" />
<p class="caption">(\#fig:sampleClustering)Figure S6: Hierarchical clustering of the samples.</p>
</div>

We can observe that samples cluster primarily by sample type, tumor or normal. TSS seems to have
a stronger effect among the normal samples, while it distributes better among the tumor samples.
We may consider discarding samples leading to an unbalanced distribution of the outcome across batches.

In Figure \@ref(fig:mdsPlot) we show the corresponding MDS plot. Here we see more clearly that the
first source of variation separates tumor from normal samples. We can also observe that two tumor
samples, corresponding to individuals `KL-8404` and `KN-8427` are separated from the rest, just as
it happens in the hierchical clustering. A closer examination of their corresponding MA-plots also
reveals a slight dependence of expression changes on average expression. We may consider discarding
these two samples and doing the MDS plot again to have a closer look to the differences among the rest
of the samples and their relationship with TSS.


```r
plotMDS(dge.filt, labels=outcome, col=batch)
legend("bottomleft", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)), inset=0.05)
```

<div class="figure" style="text-align: center">
<img src="QAanalysis_files/figure-html/mdsPlot-1.png" alt="Figure S7: Multidimensional scaling plot of the samples." width="1400" />
<p class="caption">(\#fig:mdsPlot)Figure S7: Multidimensional scaling plot of the samples.</p>
</div>

## Session information


```r
sessionInfo()
```

```
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS: /opt/R/R-3.5.0/lib64/R/lib/libRblas.so
LAPACK: /opt/R/R-3.5.0/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF8       LC_NUMERIC=C             
 [3] LC_TIME=en_US.UTF8        LC_COLLATE=en_US.UTF8    
 [5] LC_MONETARY=en_US.UTF8    LC_MESSAGES=en_US.UTF8   
 [7] LC_PAPER=en_US.UTF8       LC_NAME=C                
 [9] LC_ADDRESS=C              LC_TELEPHONE=C           
[11] LC_MEASUREMENT=en_US.UTF8 LC_IDENTIFICATION=C      

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] geneplotter_1.58.0          annotate_1.58.0            
 [3] XML_3.98-1.16               AnnotationDbi_1.42.1       
 [5] lattice_0.20-35             edgeR_3.22.5               
 [7] limma_3.36.5                SummarizedExperiment_1.10.1
 [9] DelayedArray_0.6.6          BiocParallel_1.14.2        
[11] matrixStats_0.54.0          Biobase_2.40.0             
[13] GenomicRanges_1.32.7        GenomeInfoDb_1.16.0        
[15] IRanges_2.14.12             S4Vectors_0.18.3           
[17] BiocGenerics_0.26.0         knitr_1.20                 
[19] BiocStyle_2.8.2            

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19           RColorBrewer_1.1-2     compiler_3.5.0        
 [4] highr_0.7              XVector_0.20.0         bitops_1.0-6          
 [7] tools_3.5.0            zlibbioc_1.26.0        bit_1.1-14            
[10] digest_0.6.18          memoise_1.1.0          RSQLite_2.1.1         
[13] evaluate_0.12          Matrix_1.2-14          DBI_1.0.0             
[16] yaml_2.2.0             xfun_0.3               GenomeInfoDbData_1.1.0
[19] stringr_1.3.1          bit64_0.9-7            locfit_1.5-9.1        
[22] rprojroot_1.3-2        grid_3.5.0             rmarkdown_1.10        
[25] bookdown_0.7           blob_1.1.1             magrittr_1.5          
[28] backports_1.1.2        codetools_0.2-15       htmltools_0.3.6       
[31] xtable_1.8-3           KernSmooth_2.23-15     stringi_1.2.4         
[34] RCurl_1.95-4.11       
```
