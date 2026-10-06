| cell | series | η compared / viol. | lag/T mismatch | d = ndiffs | d ≠ ndiffs at tie / defect | C-S′ worse % [W] (cases) | R step=ex % | MASE ours/R-auto med [CI] |
|---|---|---|---|---|---|---|---|---|
| wn|200 | 100 | 107 / 0 | 0 | 100 | 0 / 0 | 3.0 [1.0, 8.5] (3) | 83.0 | 1.000 [1.000, 1.000] |
| wn|500 | 100 | 103 / 0 | 0 | 100 | 0 / 0 | 1.0 [0.2, 5.4] (1) | 85.0 | 1.000 [1.000, 1.000] |
| ar1|200 | 100 | 121 / 0 | 0 | 100 | 0 / 0 | 7.0 [3.4, 13.7] (7) | 72.0 | 1.000 [1.000, 1.000] |
| ar1|500 | 100 | 117 / 0 | 0 | 100 | 0 / 0 | 12.0 [7.0, 19.8] (12) | 79.0 | 1.000 [1.000, 1.000] |
| ma1|200 | 100 | 110 / 0 | 0 | 100 | 0 / 0 | 5.0 [2.2, 11.2] (5) | 72.0 | 1.000 [1.000, 1.000] |
| ma1|500 | 100 | 109 / 0 | 0 | 100 | 0 / 0 | 5.0 [2.2, 11.2] (5) | 76.0 | 1.000 [1.000, 1.000] |
| arma11|200 | 100 | 110 / 0 | 0 | 100 | 0 / 0 | 17.0 [10.9, 25.5] (17) | 68.0 | 1.000 [1.000, 1.000] |
| arma11|500 | 100 | 109 / 0 | 0 | 100 | 0 / 0 | 15.0 [9.3, 23.3] (15) | 63.0 | 1.000 [1.000, 1.000] |
| ar2|200 | 100 | 139 / 0 | 0 | 100 | 0 / 0 | 26.3 [18.6, 35.7] (26) | 54.0 | 1.000 [1.000, 1.000] |
| ar2|500 | 100 | 130 / 0 | 0 | 100 | 0 / 0 | 14.0 [8.5, 22.1] (14) | 44.0 | 1.000 [1.000, 1.000] |
| ar098|200 | 100 | 192 / 0 | 0 | 100 | 0 / 0 | 0.0 [0.0, 3.8] (0) | 65.0 | 1.000 [1.000, 1.000] |
| ar098|500 | 100 | 195 / 0 | 0 | 100 | 0 / 0 | 2.0 [0.6, 7.1] (2) | 70.0 | 1.000 [1.000, 1.000] |
| ima|200 | 100 | 199 / 0 | 0 | 100 | 0 / 0 | 3.0 [1.0, 8.5] (3) | 64.0 | 1.000 [1.000, 1.000] |
| ima|500 | 100 | 209 / 0 | 0 | 100 | 0 / 0 | 8.0 [4.1, 15.0] (8) | 64.0 | 1.000 [1.000, 1.000] |
| ari|200 | 100 | 206 / 0 | 0 | 100 | 0 / 0 | 6.0 [2.8, 12.5] (6) | 75.0 | 1.000 [1.000, 1.000] |
| ari|500 | 100 | 203 / 0 | 0 | 100 | 0 / 0 | 6.1 [2.8, 12.6] (6) | 78.0 | 1.000 [1.000, 1.000] |
| arima111|200 | 100 | 208 / 0 | 0 | 100 | 0 / 0 | 28.0 [20.1, 37.5] (28) | 72.0 | 1.000 [1.000, 1.000] |
| arima111|500 | 100 | 209 / 0 | 0 | 100 | 0 / 0 | 21.0 [14.2, 30.0] (21) | 67.0 | 1.000 [1.000, 1.000] |
| imadrift|200 | 100 | 204 / 0 | 0 | 100 | 0 / 0 | 3.0 [1.0, 8.5] (3) | 84.0 | 1.000 [1.000, 1.000] |
| imadrift|500 | 100 | 206 / 0 | 0 | 100 | 0 / 0 | 3.0 [1.0, 8.5] (3) | 66.0 | 1.000 [1.000, 1.000] |
| airline|200 | 100 | 179 / 0 | 0 | 100 | 0 / 0 | 12.0 [7.0, 19.8] (12) | 73.0 | 1.000 [1.000, 1.000] |
| airline|500 | 100 | 198 / 0 | 0 | 100 | 0 / 0 | 7.0 [3.4, 13.7] (7) | 80.0 | 1.000 [1.000, 1.000] |

F_S classical vs STL (seasonal-period series):
D* = nsdiffs_R on 193 / 200

DEFECTS (C-d′ η/lag/d): none
STATISTICAL GATES (C-S′, MASE vs R): FAIL

details:
C-S′: wn_n200_r6 ours 566.8290 (2,2,0,0,c=1) vs R order (c: true, p: 4, q: 0, sp: 0, sq: 0) 565.5318
C-S′: wn_n200_r25 ours 583.8167 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 583.1679
C-S′: wn_n200_r47 ours 582.9473 (1,2,0,0,c=1) vs R order (c: true, p: 3, q: 3, sp: 0, sq: 0) 579.0703
C-S′: wn_n500_r27 ours 1414.6214 (3,1,0,0,c=1) vs R order (c: true, p: 1, q: 3, sp: 0, sq: 0) 1413.4762
C-S′: ar1_n200_r1 ours 574.0995 (0,4,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 565.0752
C-S′: ar1_n200_r2 ours 594.6711 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 594.6394
C-S′: ar1_n200_r20 ours 563.8615 (0,3,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 561.9648
C-S′: ar1_n200_r25 ours 560.4861 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 560.4611
C-S′: ar1_n200_r52 ours 573.3448 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 572.9218
C-S′: ar1_n200_r57 ours 589.3411 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 589.1702
C-S′: ar1_n200_r96 ours 547.1956 (0,3,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 543.5807
C-S′: ar1_n500_r9 ours 1363.4529 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1363.4074
C-S′: ar1_n500_r10 ours 1485.3593 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1485.1783
C-S′: ar1_n500_r32 ours 1408.3259 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1407.8717
C-S′: ar1_n500_r50 ours 1478.2403 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1478.1063
C-S′: ar1_n500_r54 ours 1458.6235 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1458.5816
C-S′: ar1_n500_r60 ours 1456.1837 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1456.0268
C-S′: ar1_n500_r61 ours 1408.4611 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1408.0414
C-S′: ar1_n500_r74 ours 1442.3062 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1442.2550
C-S′: ar1_n500_r76 ours 1387.0456 (2,2,0,0,c=1) vs R order (c: true, p: 1, q: 3, sp: 0, sq: 0) 1385.9546
C-S′: ar1_n500_r92 ours 1399.3373 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 0) 1398.8481
C-S′: ar1_n500_r95 ours 1446.9646 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1446.8459
C-S′: ar1_n500_r99 ours 1458.5743 (1,4,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1457.9741
C-S′: ma1_n200_r36 ours 596.3976 (3,1,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 594.2524
C-S′: ma1_n200_r40 ours 548.2608 (2,2,0,0,c=1) vs R order (c: true, p: 0, q: 4, sp: 0, sq: 0) 547.2956
C-S′: ma1_n200_r55 ours 600.3582 (2,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 596.9385
C-S′: ma1_n200_r85 ours 556.5275 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 556.1206
C-S′: ma1_n200_r86 ours 592.5910 (2,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 589.9261
C-S′: ma1_n500_r16 ours 1431.8263 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 1431.2559
C-S′: ma1_n500_r43 ours 1419.0409 (2,2,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 1418.3652
C-S′: ma1_n500_r67 ours 1440.4627 (1,2,0,0,c=1) vs R order (c: true, p: 0, q: 3, sp: 0, sq: 0) 1440.2154
C-S′: ma1_n500_r74 ours 1392.3189 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 1392.0442
C-S′: ma1_n500_r88 ours 1401.4459 (4,1,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 1399.9714
C-S′: arma11_n200_r13 ours 572.6411 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 572.4306
C-S′: arma11_n200_r25 ours 562.6078 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 560.0690
C-S′: arma11_n200_r26 ours 575.1450 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 574.1313
C-S′: arma11_n200_r31 ours 549.9845 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 549.4810
C-S′: arma11_n200_r32 ours 573.7513 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 572.4965
C-S′: arma11_n200_r33 ours 568.2376 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 566.9309
C-S′: arma11_n200_r35 ours 556.3103 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 554.2447
C-S′: arma11_n200_r46 ours 575.8936 (2,1,0,0,c=0) vs R order (c: false, p: 1, q: 2, sp: 0, sq: 0) 575.2765
C-S′: arma11_n200_r53 ours 575.1317 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 574.8046
C-S′: arma11_n200_r73 ours 550.2397 (2,2,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 550.2027
C-S′: arma11_n200_r76 ours 571.2075 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 571.1941
C-S′: arma11_n200_r77 ours 586.0003 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 585.0857
C-S′: arma11_n200_r78 ours 568.9730 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 566.8470
C-S′: arma11_n200_r85 ours 553.0713 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 552.2668
C-S′: arma11_n200_r86 ours 591.6760 (1,2,0,0,c=1) vs R order (c: true, p: 0, q: 3, sp: 0, sq: 0) 590.9259
C-S′: arma11_n200_r91 ours 592.0354 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 591.7608
C-S′: arma11_n200_r97 ours 576.7448 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 575.8301
C-S′: arma11_n500_r15 ours 1467.8017 (3,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1465.0574
C-S′: arma11_n500_r17 ours 1493.1448 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1493.1087
C-S′: arma11_n500_r18 ours 1477.7062 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1476.4952
C-S′: arma11_n500_r20 ours 1406.2121 (3,1,0,0,c=0) vs R order (c: false, p: 1, q: 3, sp: 0, sq: 0) 1406.0910
C-S′: arma11_n500_r29 ours 1440.0367 (3,1,0,0,c=0) vs R order (c: false, p: 1, q: 2, sp: 0, sq: 0) 1438.5412
C-S′: arma11_n500_r38 ours 1428.8780 (1,2,0,0,c=1) vs R order (c: true, p: 0, q: 3, sp: 0, sq: 0) 1428.6168
C-S′: arma11_n500_r40 ours 1448.7355 (2,2,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 1448.3368
C-S′: arma11_n500_r43 ours 1396.7516 (1,2,0,0,c=1) vs R order (c: true, p: 0, q: 3, sp: 0, sq: 0) 1395.1537
C-S′: arma11_n500_r44 ours 1386.9993 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1386.9926
C-S′: arma11_n500_r58 ours 1453.4927 (3,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1451.2777
C-S′: arma11_n500_r65 ours 1424.3186 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 1, sp: 0, sq: 0) 1423.8963
C-S′: arma11_n500_r77 ours 1400.9273 (4,1,0,0,c=0) vs R order (c: false, p: 0, q: 4, sp: 0, sq: 0) 1400.6946
C-S′: arma11_n500_r91 ours 1437.2563 (3,1,0,0,c=0) vs R order (c: false, p: 2, q: 3, sp: 0, sq: 0) 1434.7703
C-S′: arma11_n500_r98 ours 1407.5550 (2,2,0,0,c=1) vs R order (c: true, p: 2, q: 4, sp: 0, sq: 0) 1406.2386
C-S′: arma11_n500_r99 ours 1400.5126 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 1, sp: 0, sq: 0) 1400.0510
C-S′: ar2_n200_r0 ours 564.7362 (2,2,0,0,c=1) vs R order (c: true, p: 1, q: 3, sp: 0, sq: 0) 563.1398
C-S′: ar2_n200_r2 ours 582.5393 (1,1,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 581.5329
C-S′: ar2_n200_r3 ours 597.1693 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 596.4185
C-S′: ar2_n200_r9 ours 570.0873 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 569.9880
C-S′: ar2_n200_r20 ours 592.6196 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 592.3592
C-S′: ar2_n200_r23 ours 564.7475 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 563.8247
C-S′ skipped: ar2_n200_r24 R step order (c: true, p: 2, q: 2, sp: 0, sq: 0) rejectedRootsNearUnit with us
C-S′: ar2_n200_r26 ours 574.2207 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 572.4425
C-S′: ar2_n200_r35 ours 568.3906 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 568.3214
C-S′: ar2_n200_r36 ours 559.6411 (1,2,0,0,c=1) vs R order (c: true, p: 3, q: 0, sp: 0, sq: 0) 558.6533
C-S′: ar2_n200_r41 ours 566.2069 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 565.1363
C-S′: ar2_n200_r42 ours 572.8434 (2,1,0,0,c=0) vs R order (c: false, p: 1, q: 2, sp: 0, sq: 0) 572.4207
C-S′: ar2_n200_r43 ours 546.4043 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 544.7962
C-S′: ar2_n200_r49 ours 550.9523 (2,2,0,0,c=1) vs R order (c: true, p: 1, q: 3, sp: 0, sq: 0) 550.3264
C-S′: ar2_n200_r56 ours 637.6289 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 635.0702
C-S′: ar2_n200_r58 ours 564.2737 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 561.8067
C-S′: ar2_n200_r62 ours 551.1963 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 550.3089
C-S′: ar2_n200_r69 ours 550.4076 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 548.0635
C-S′: ar2_n200_r74 ours 583.6980 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 583.5929
C-S′: ar2_n200_r76 ours 592.8795 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 590.1544
C-S′: ar2_n200_r78 ours 543.8479 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 0) 543.3829
C-S′: ar2_n200_r82 ours 576.2315 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 575.4240
C-S′: ar2_n200_r89 ours 563.0394 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 562.1689
C-S′: ar2_n200_r93 ours 546.3382 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 544.9725
C-S′: ar2_n200_r94 ours 548.4851 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 547.6730
C-S′: ar2_n200_r97 ours 561.1970 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 560.6349
C-S′: ar2_n200_r99 ours 580.6622 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 578.6788
C-S′: ar2_n500_r4 ours 1418.4507 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1417.1965
C-S′: ar2_n500_r12 ours 1371.7414 (2,1,0,0,c=1) vs R order (c: true, p: 3, q: 0, sp: 0, sq: 0) 1371.2952
C-S′: ar2_n500_r19 ours 1393.4089 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1392.5671
C-S′: ar2_n500_r20 ours 1469.3384 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1467.1919
C-S′: ar2_n500_r47 ours 1457.7858 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1454.9170
C-S′: ar2_n500_r67 ours 1385.5380 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1383.6735
C-S′: ar2_n500_r76 ours 1426.9077 (2,1,0,0,c=1) vs R order (c: true, p: 3, q: 0, sp: 0, sq: 0) 1426.4582
C-S′: ar2_n500_r77 ours 1463.7146 (1,2,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 1459.2692
C-S′: ar2_n500_r82 ours 1425.7528 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1424.6985
C-S′: ar2_n500_r87 ours 1467.4470 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1465.0629
C-S′: ar2_n500_r91 ours 1413.7254 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1413.1441
C-S′: ar2_n500_r97 ours 1481.8083 (1,2,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 1476.1048
C-S′: ar2_n500_r98 ours 1426.2309 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 1424.6844
C-S′: ar2_n500_r99 ours 1441.9780 (1,1,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 1437.3266
C-S′ skipped: ar098_n200_r59 R step order (c: true, p: 1, q: 1, sp: 0, sq: 0) rejectedRootsNearUnit with us
C-S′ skipped: ar098_n200_r65 R step order (c: false, p: 1, q: 2, sp: 0, sq: 0) rejectedRootsNearUnit with us
C-S′: ar098_n500_r75 ours 1478.5075 (2,2,0,0,c=1) vs R order (c: true, p: 3, q: 1, sp: 0, sq: 0) 1478.4920
C-S′ skipped: ar098_n500_r77 R step order (c: false, p: 2, q: 2, sp: 0, sq: 0) rejectedRootsNearUnit with us
C-S′: ar098_n500_r95 ours 1395.6847 (1,1,0,0,c=0) vs R order (c: false, p: 3, q: 1, sp: 0, sq: 0) 1394.4999
C-S′: ima_n200_r25 ours 581.1654 (2,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 580.1569
C-S′: ima_n200_r40 ours 561.8848 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 561.2118
C-S′: ima_n200_r55 ours 571.9067 (1,1,0,0,c=1) vs R order (c: false, p: 3, q: 0, sp: 0, sq: 0) 564.7228
C-S′: ima_n500_r39 ours 1352.0351 (1,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 1351.7437
C-S′: ima_n500_r60 ours 1377.0369 (1,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 1376.9502
C-S′: ima_n500_r67 ours 1452.2436 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 1452.1064
C-S′: ima_n500_r71 ours 1417.0051 (1,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 1416.6551
C-S′: ima_n500_r77 ours 1445.7438 (1,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 1445.7168
C-S′: ima_n500_r85 ours 1369.9154 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 0) 1369.7835
C-S′: ima_n500_r98 ours 1379.2382 (2,1,0,0,c=0) vs R order (c: false, p: 3, q: 0, sp: 0, sq: 0) 1378.7022
C-S′: ima_n500_r99 ours 1427.1406 (1,1,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 1426.5363
C-S′: ari_n200_r7 ours 549.1648 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 549.1227
C-S′: ari_n200_r43 ours 580.3698 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 579.8740
C-S′: ari_n200_r56 ours 621.5228 (2,0,0,0,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 0) 621.2068
C-S′: ari_n200_r67 ours 528.6683 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 528.3560
C-S′: ari_n200_r75 ours 577.8952 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 577.7714
C-S′: ari_n200_r87 ours 559.4075 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 558.1682
C-S′: ari_n500_r1 ours 1411.2187 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1410.9016
C-S′: ari_n500_r6 ours 1447.8066 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1447.6521
C-S′: ari_n500_r25 ours 1374.5070 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1374.2291
C-S′: ari_n500_r35 ours 1362.8654 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1362.8393
C-S′: ari_n500_r68 ours 1389.6667 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1389.6644
C-S′: ari_n500_r81 ours 1435.5276 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1435.4885
C-S′ skipped: ari_n500_r89 R step order (c: true, p: 2, q: 1, sp: 0, sq: 0) rejectedRootsNearUnit with us
C-S′: arima111_n200_r3 ours 532.0706 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 531.0077
C-S′: arima111_n200_r4 ours 563.8606 (2,2,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 563.2642
C-S′: arima111_n200_r7 ours 552.7734 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 551.5601
C-S′: arima111_n200_r10 ours 559.6401 (0,2,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 558.9085
C-S′: arima111_n200_r12 ours 569.7047 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 568.6437
C-S′: arima111_n200_r13 ours 528.5678 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 527.0497
C-S′: arima111_n200_r16 ours 599.0285 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 598.6500
C-S′: arima111_n200_r26 ours 569.3162 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 568.7192
C-S′: arima111_n200_r27 ours 551.2789 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 550.0245
C-S′: arima111_n200_r30 ours 531.8182 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 0) 531.0821
C-S′: arima111_n200_r31 ours 548.3654 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 547.4309
C-S′: arima111_n200_r32 ours 550.3911 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 549.9398
C-S′: arima111_n200_r33 ours 570.5461 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 569.5140
C-S′: arima111_n200_r36 ours 550.5600 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 550.4143
C-S′: arima111_n200_r39 ours 569.4620 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 568.7643
C-S′: arima111_n200_r41 ours 583.5164 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 582.4340
C-S′: arima111_n200_r51 ours 552.0864 (1,1,0,0,c=0) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 550.8595
C-S′: arima111_n200_r54 ours 563.2715 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 563.0516
C-S′: arima111_n200_r57 ours 570.9914 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 570.4275
C-S′: arima111_n200_r58 ours 553.4449 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 551.6869
C-S′: arima111_n200_r59 ours 552.5591 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 551.9783
C-S′: arima111_n200_r60 ours 564.0937 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 563.8686
C-S′: arima111_n200_r61 ours 540.9050 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 540.7374
C-S′: arima111_n200_r67 ours 561.0478 (0,2,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 559.8155
C-S′: arima111_n200_r73 ours 538.0260 (1,1,0,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 0) 537.5826
C-S′: arima111_n200_r75 ours 574.0233 (2,3,0,0,c=0) vs R order (c: false, p: 2, q: 5, sp: 0, sq: 0) 571.9014
C-S′: arima111_n200_r89 ours 600.6557 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 600.3730
C-S′: arima111_n200_r99 ours 578.7388 (1,1,0,0,c=1) vs R order (c: true, p: 4, q: 1, sp: 0, sq: 0) 575.8287
C-S′: arima111_n500_r7 ours 1460.7580 (2,2,0,0,c=0) vs R order (c: false, p: 3, q: 1, sp: 0, sq: 0) 1460.6869
C-S′: arima111_n500_r9 ours 1434.7229 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1434.4881
C-S′: arima111_n500_r16 ours 1408.7880 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1408.2931
C-S′: arima111_n500_r17 ours 1410.4261 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1409.3582
C-S′: arima111_n500_r28 ours 1460.0296 (2,0,0,0,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 0) 1459.5353
C-S′: arima111_n500_r31 ours 1438.2199 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1437.9964
C-S′: arima111_n500_r33 ours 1377.0941 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1376.7620
C-S′: arima111_n500_r36 ours 1409.0219 (3,1,0,0,c=0) vs R order (c: false, p: 2, q: 2, sp: 0, sq: 0) 1408.6129
C-S′: arima111_n500_r38 ours 1412.0238 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1411.3273
C-S′: arima111_n500_r49 ours 1318.3702 (2,0,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1317.1238
C-S′: arima111_n500_r52 ours 1427.0034 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1425.9683
C-S′: arima111_n500_r55 ours 1405.1438 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 1403.5661
C-S′: arima111_n500_r56 ours 1379.6809 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 0) 1377.9209
C-S′: arima111_n500_r58 ours 1414.6209 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1413.7153
C-S′: arima111_n500_r64 ours 1386.5881 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 1385.3235
C-S′: arima111_n500_r72 ours 1394.6214 (1,1,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1393.8874
C-S′: arima111_n500_r81 ours 1411.5080 (0,2,0,0,c=1) vs R order (c: true, p: 1, q: 1, sp: 0, sq: 0) 1410.5610
C-S′: arima111_n500_r83 ours 1386.1591 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 0, sp: 0, sq: 0) 1383.3562
C-S′: arima111_n500_r93 ours 1406.0678 (1,2,0,0,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 0) 1405.5780
C-S′: arima111_n500_r96 ours 1354.9557 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 1354.9480
C-S′: arima111_n500_r99 ours 1415.0155 (1,2,0,0,c=0) vs R order (c: false, p: 0, q: 3, sp: 0, sq: 0) 1413.9707
C-S′: imadrift_n200_r12 ours 549.6687 (1,2,0,0,c=1) vs R order (c: true, p: 2, q: 1, sp: 0, sq: 0) 549.1038
C-S′: imadrift_n200_r45 ours 553.0416 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 552.9203
C-S′: imadrift_n200_r61 ours 578.5514 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 578.0754
C-S′: imadrift_n500_r40 ours 1439.1306 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 1438.5386
C-S′: imadrift_n500_r59 ours 1454.8267 (2,3,0,0,c=1) vs R order (c: true, p: 1, q: 4, sp: 0, sq: 0) 1454.4217
C-S′: imadrift_n500_r60 ours 1420.9244 (1,1,0,0,c=1) vs R order (c: true, p: 0, q: 2, sp: 0, sq: 0) 1420.4605
C-S′: airline_n200_r6 ours 528.5535 (1,2,0,1,c=0) vs R order (c: false, p: 2, q: 1, sp: 0, sq: 1) 526.9512
C-S′: airline_n200_r18 ours 581.9256 (1,1,0,1,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 1) 581.2878
C-S′: airline_n200_r21 ours 524.7087 (2,3,1,1,c=0) vs R order (c: false, p: 2, q: 3, sp: 0, sq: 2) 523.9202
C-S′: airline_n200_r22 ours 538.0431 (2,0,2,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 1) 535.2035
C-S′: airline_n200_r49 ours 561.8771 (1,1,1,1,c=0) vs R order (c: false, p: 1, q: 1, sp: 0, sq: 2) 561.1020
C-S′: airline_n200_r59 ours 517.0803 (1,1,0,1,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 1) 515.7066
C-S′: airline_n200_r69 ours 513.7959 (1,1,0,1,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 1) 513.1044
C-S′: airline_n200_r74 ours 553.5931 (0,2,1,1,c=0) vs R order (c: false, p: 1, q: 1, sp: 1, sq: 1) 551.4997
C-S′: airline_n200_r78 ours 533.5460 (2,0,1,0,c=1) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 1) 529.0813
C-S′: airline_n200_r83 ours 542.8063 (1,1,0,1,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 1) 542.0282
C-S′: airline_n200_r86 ours 513.8638 (1,1,1,1,c=1) vs R order (c: true, p: 3, q: 2, sp: 0, sq: 2) 512.3787
C-S′: airline_n200_r99 ours 530.6557 (1,1,0,1,c=0) vs R order (c: true, p: 2, q: 0, sp: 0, sq: 1) 529.9722
C-S′: airline_n500_r16 ours 1393.5169 (0,1,0,1,c=0) vs R order (c: false, p: 1, q: 0, sp: 0, sq: 1) 1392.7090
C-S′: airline_n500_r17 ours 1371.8205 (1,1,0,2,c=0) vs R order (c: false, p: 0, q: 2, sp: 0, sq: 2) 1371.3127
C-S′: airline_n500_r21 ours 1416.3780 (0,1,1,1,c=0) vs R order (c: false, p: 0, q: 1, sp: 0, sq: 2) 1416.1849
C-S′: airline_n500_r34 ours 1546.8080 (3,2,0,2,c=0) vs R order (c: false, p: 5, q: 0, sp: 2, sq: 0) 1513.5736
C-S′: airline_n500_r42 ours 1336.2674 (2,3,1,1,c=0) vs R order (c: false, p: 1, q: 2, sp: 0, sq: 1) 1334.4568
C-S′: airline_n500_r52 ours 1365.6277 (0,1,1,1,c=0) vs R order (c: false, p: 0, q: 1, sp: 0, sq: 2) 1365.4664
C-S′: airline_n500_r60 ours 1384.6309 (0,1,1,1,c=0) vs R order (c: false, p: 0, q: 1, sp: 0, sq: 2) 1384.3888