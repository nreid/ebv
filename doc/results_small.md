# Results

## Small RNA libraries

I analyzed the data using two pieces of software:

`ShortStack`

    Axtell MJ. (2013) ShortStack: Comprehensive annotation and
    quantification of small RNA genes. RNA 19:740-751.
    doi:10.1261/rna.035279.112

`miRDeep2`

    Friedländer, Marc R., et al. "miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades." Nucleic acids research 40.1 (2012): 37-52.

These are very different approaches. `ShortStack` was originally developed for plant data, but is not limited to analyzing only miRNA data. It is included here because it is simple to use and well documented. It did seem to have some flaws, however. At least one pair of what appear to be clearly different miRNA pileups have been subsumed into a single locus. I'm guessing it is not great at identifying distinct loci whose pileups may partially overlap when sequencing is extremely deep, as it is here. 

`miRDeep2` is a bit more opaque, difficult to run, and the output files were harder to parse. It also only discoveres and quantifies miRNA, which seem to be a relatively small proportion of the total small RNA in these libraries. It takes in known a priori miRNAs as well as discovers new ones. 


| Sample            |    total | shortstack | mirdeep2 |shortstack_pct | mirdeep2_pct |
| ------------------|----------|------------|----------|---------------|------------- |
| Daudi_Ago1KD_S2   | 37690751 |   36684233 |    86272 |     0.9732954 | 0.0022889435 |
| Daudi_Ago2KD_S8   | 37838612 |   37139800 |    12498 |     0.9815318 | 0.0003302975 |
| Daudi_Ago3KD_S1   | 35847271 |   35061107 |    52443 |     0.9780691 | 0.0014629566 |
| Daudi_Ago4KD_S7   | 38849754 |   38064144 |    93806 |     0.9797783 | 0.0024145841 |
| Daudi_CONTROL_S3  | 35911357 |   35296735 |    81018 |     0.9828850 | 0.0022560551 |
| Daudi_DicerKD_S4  | 31137968 |   30190821 |   216044 |     0.9695822 | 0.0069382819 |
| Daudi_DroshaKD_S5 | 42491783 |   41406032 |    65306 |     0.9744480 | 0.0015369089 |
| Daudi_LaKD_S6     | 38041319 |   37249392 |   139916 |     0.9791825 | 0.0036780008 |
| Hsa_Ago1KD_S6     | 45904528 |   44096436 |    33107 |     0.9606119 | 0.0007212143 |
| Hsa_Ago2KD_S8     | 34400525 |   33543241 |    40541 |     0.9750793 | 0.0011784995 |
| Hsa_Ago3KD_S5     | 43919512 |   42003603 |   138173 |     0.9563768 | 0.0031460504 |
| Hsa_Ago4KD_S3     | 32877686 |   32036686 |    32369 |     0.9744203 | 0.0009845279 |
| Hsa_CONTROL_S1    | 34269714 |   33214043 |   290339 |     0.9691952 | 0.0084721746 |
| Hsa_DicerKD_S10   | 31011440 |   30164662 |    22851 |     0.9726947 | 0.0007368571 |
| Hsa_DroshaKD_S9   | 43649632 |   41153508 |    58077 |     0.9428145 | 0.0013305267 |
| Hsa_LaKD_S7       | 38900536 |   37884251 |    12791 |     0.9738748 | 0.0003288130 |
| RE_Ago1KD_S4      | 12024805 |   11788186 |   204583 |     0.9803224 | 0.0170134152 |
| RE_Ago2KD_S5      | 21599261 |   21287827 |   102015 |     0.9855813 | 0.0047230783 |
| RE_Ago3KD_S8      | 20532202 |   20121683 |    35757 |     0.9800061 | 0.0017415083 |
| RE_Ago4KD_S7      | 19044265 |   18704282 |    44222 |     0.9821477 | 0.0023220639 |
| RE_CONTROL_S7     | 24530322 |   23386737 |   102144 |     0.9533808 | 0.0041639894 |
| RE_DicerKD_S6     | 19794234 |   19223768 |    92230 |     0.9711802 | 0.0046594377 |
| RE_DroshaKD_S1    | 20013765 |   19420564 |   116512 |     0.9703603 | 0.0058215933 |
| RE_LaKD_S3        | 16808760 |   16308042 |   150788 |     0.9702109 | 0.0089707986 |
| SNU_ago1KD_S5     | 57290875 |   54765836 |  4303371 |     0.9559260 | 0.0751144227 |
| SNU_ago2KD_S6     | 55376358 |   53042587 |  1008138 |     0.9578562 | 0.0182052059 |
| SNU_ago3KD_S7     | 39212773 |   37398250 |  1261676 |     0.9537262 | 0.0321751282 |
| SNU_ago4KD_S4     | 38611338 |   36904035 |  2216997 |     0.9557823 | 0.0574182899 |
| SNU_CONTROL_X     | 23358020 |   22907909 |   342292 |     0.9807299 | 0.0146541530 |
| SNU_dicerKD_S8    | 55815305 |   53622247 |  2535480 |     0.9607087 | 0.0454262500 |
| SNU_droshaKD_S3   | 63704620 |   61334249 |   410972 |     0.9627912 | 0.0064512119 |
| SNU_laKD_S1       | 40057104 |   38326299 |  1367237 |     0.9567916 | 0.0341321979 |



### Total reads analyzed, EBV only

| Sample            |    total | shortstack | mirdeep2 | shortstack_% | mirdeep2_%|
| ------------------|----------|------------|----------|----------------|-------------|
| Daudi_Ago1KD_S2   | 37690751 |      30835 |     2220 |   8.181052e-04 | 5.890039e-05|
| Daudi_Ago2KD_S8   | 37838612 |       9691 |      359 |   2.561140e-04 | 9.487663e-06|
| Daudi_Ago3KD_S1   | 35847271 |      20565 |     1572 |   5.736838e-04 | 4.385271e-05|
| Daudi_Ago4KD_S7   | 38849754 |      27455 |     2843 |   7.066969e-04 | 7.317936e-05|
| Daudi_CONTROL_S3  | 35911357 |      22619 |     4017 |   6.298565e-04 | 1.118588e-04|
| Daudi_DicerKD_S4  | 31137968 |      53111 |     6483 |   1.705667e-03 | 2.082024e-04|
| Daudi_DroshaKD_S5 | 42491783 |      29287 |     1592 |   6.892391e-04 | 3.746607e-05|
| Daudi_LaKD_S6     | 38041319 |      37182 |     4027 |   9.774109e-04 | 1.058586e-04|
| Hsa_Ago1KD_S6     | 45904528 |      92429 |     1012 |   2.013505e-03 | 2.204576e-05|
| Hsa_Ago2KD_S8     | 34400525 |     132276 |     1628 |   3.845174e-03 | 4.732486e-05|
| Hsa_Ago3KD_S5     | 43919512 |     114078 |     2764 |   2.597433e-03 | 6.293330e-05|
| Hsa_Ago4KD_S3     | 32877686 |      96807 |     1360 |   2.944459e-03 | 4.136544e-05|
| Hsa_CONTROL_S1    | 34269714 |      86390 |     7931 |   2.520885e-03 | 2.314288e-04|
| Hsa_DicerKD_S10   | 31011440 |      16162 |      595 |   5.211625e-04 | 1.918647e-05|
| Hsa_DroshaKD_S9   | 43649632 |      94038 |     1884 |   2.154382e-03 | 4.316188e-05|
| Hsa_LaKD_S7       | 38900536 |     101071 |      350 |   2.598190e-03 | 8.997305e-06|
| RE_Ago1KD_S4      | 12024805 |         29 |        1 |   2.411682e-06 | 8.316143e-08|
| RE_Ago2KD_S5      | 21599261 |          0 |        0 |   0.000000e+00 | 0.000000e+00|
| RE_Ago3KD_S8      | 20532202 |          1 |        0 |   4.870398e-08 | 0.000000e+00|
| RE_Ago4KD_S7      | 19044265 |          0 |        0 |   0.000000e+00 | 0.000000e+00|
| RE_CONTROL_S7     | 24530322 |          3 |        0 |   1.222976e-07 | 0.000000e+00|
| RE_DicerKD_S6     | 19794234 |          0 |        0 |   0.000000e+00 | 0.000000e+00|
| RE_DroshaKD_S1    | 20013765 |        446 |       22 |   2.228466e-05 | 1.099243e-06|
| RE_LaKD_S3        | 16808760 |          4 |        0 |   2.379712e-07 | 0.000000e+00|
| SNU_ago1KD_S5     | 57290875 |    1518108 |   505480 |   2.649825e-02 | 8.823046e-03|
| SNU_ago2KD_S6     | 55376358 |     432637 |   129816 |   7.812666e-03 | 2.344250e-03|
| SNU_ago3KD_S7     | 39212773 |     551122 |   155093 |   1.405466e-02 | 3.955165e-03|
| SNU_ago4KD_S4     | 38611338 |     733124 |   244157 |   1.898727e-02 | 6.323453e-03|
| SNU_CONTROL_X     | 23358020 |     108273 |    30195 |   4.635367e-03 | 1.292704e-03|
| SNU_dicerKD_S8    | 55815305 |     861734 |   275588 |   1.543903e-02 | 4.937499e-03|
| SNU_droshaKD_S3   | 63704620 |     228712 |    46797 |   3.590195e-03 | 7.345935e-04|
| SNU_laKD_S1       | 40057104 |     402880 |   121341 |   1.005764e-02 | 3.029201e-03|