#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Y-chromosome Haplogroup Classifier v1.0
Optimized for low-density SNP data with PLINK compatibility

v1.0Core improvements:
1. Combines v2.4 ancestor chain validation + v3.4 clear framework
2. Strict classification logic: traverse tree downward, validate each node
3. Supports downstream inference: downstream derived + upstream missing → can infer upstream
4. Conflict detection: downstream derived + upstream ancestral → mark warning
5. Detailed diagnostic output for review and debugging

Classification logic (simulates manual classification):
1. Start from root node, check defining SNPs at each level
2. Derived → enter node, continue downward
3. Ancestral → stop, do not enter branch
4. Missing → check downstream for derived evidence (can infer)
5. Final result = deepest validated node
"""

import os
import sys
import argparse
import gzip
from collections import defaultdict
from datetime import datetime

VERSION = "1.0"

# ============================================================
# Phylogenetic tree definition (based on ISOGG 2019-2020)
# Using indented text format for easy maintenance and review
# ============================================================
OFFICIAL_TREE = """
Y	Root
	Y	Root (Y-Adam)
		A0000	A8864
		A000-T	A8835
			A000	A10805
				A000a	A21565
				A000b	A10801
					A000b1	A10765
			A00-T	PR2921
				A00	AF6
					A00a	L1149
					A00b	A4987
					A00c	A12222
				A0-T	L1085
					A0	CTS2809.1
						A0a	L979
							A0a1	L1070
								A0a1a	P114
								A0a1b	L1289
							A0a2	L981
						A0b	L92.2
					A1	P305
						A1a	M31
						A1b	P108
							A1b1	L419
								A1b1a	L602
									A1b1a1	M14
										A1b1a1a	M6
											A1b1a1a1	P28
											A1b1a1a2	L963
												A1b1a1a2a	M114
														A1b1a1a2a1a	P262
								A1b1b	M32
									A1b1b1	M28
									A1b1b2	L427
										A1b1b2a	M51
											A1b1b2a1	P291
												A1b1b2a1a	P102
										A1b1b2b	M13
											A1b1b2b1	M118
							BT	M91
						BT	M91
							B	M60
								B1	M236
									B1a	M146
									B1b	V1108
								B2	M182
									B2a	M150
												B2a1a1	M218
														B2a1a1a1	M109
																B2a1a1a1a2	CTS4624
														B2a1a1a2	G1
												B2a1a2	M108.1
													B2a1a2a	M43
									B2b	M112
										B2b1	M192
													B2b1a1b	M30
														B2b1a1b1	M108.2
														B2b1a1c1	M211
													B2b1b1a	P6
										B2b2	P112
								B3	L1388
									B3a	Z11577
							CT	M168
								DE	M145
									D	CTS3946
										D1	M174
											D1a	CTS11577
												D1a1	F6251
													D1a1a	M15
														D1a1a1	F849
															D1a1a1a	N1
																D1a1a1a1	Z27269
																	D1a1a1a1a	PH4979
																		D1a1a1a1a1	BY15199
																			D1a1a1a1a1b	F729
																				D1a1a1a1a1b2	Y157489
																	D1a1a1a1b	Z31591
																		D1a1a1a1b1	Z31599
																			D1a1a1a1b1b	Z31584
															D1a1a1b	F1070
																D1a1a1b1	BY12793
																	D1a1a1b1a	BY12738
													D1a1b	P99
															D1a1b1a	M533
																D1a1b1a1	PH116
																	D1a1b1a1a	F1771
																		D1a1b1a1a1	BY12625
												D1a2	Z3660
													D1a2a	M64.1
														D1a2a1	M116.1
															D1a2a1a	M125
																	D1a2a1a1a	P12_1
																D1a2a1a2	IMS-JST022457
																	D1a2a1b2b	IMS-JST006841
																		D1a2a1b2b1	CTS3397
																			D1a2a1b2b1a	Z1500
																				D1a2a1b2b1a1	Z1504
																					D1a2a1b2b1a1a	CTS8093
																						D1a2a1b2b1a1a1	FGC6373
																							D1a2a1b2b1a1a1a	L137.3
																							D1a2a1b2b1a1a1b	FGC6372
																						D1a2a1b2b1a1a3	BY45234
																						D1a2a1b2b1a1a9	Z40687
																					D1a2a1b2b1a1b	BY149852
																D1a2a1a3	CTS10972
																	D1a2a1a3a	Z18556. Z18566
															D1a2a1c	CTS6609
																D1a2a1c1	CTS1897
																	D1a2a1c1a	CTS11032
																		D1a2a1c1a1	CTS218
																			D1a2a1c1a1a	CTS6909
																			D1a2a1c1a1b	CTS3033
																	D1a2a1c1b	CTS1964
																		D1a2a1c1b1	BY169023
																	D1a2a1c1c	Z30644
																D1a2a1c2	CTS103
														D1a2a2	CTS131
															D1a2a2a	CTS220
																D1a2a2a1	CTS10495
																D1a2a2a2	CTS11285
																	D1a2a2a2a	PH2316
																		D1a2a2a2a1	Z38287
																		D1a2a2a2a2	Z38289
																	D1a2a2a2b	CTS288
													D1a2b	Y34637
														D1a2b1	Y34659
															D1a2b1a	Y35020
														D1a2b2	Y34953
											D1b	L1378
									E	M96
										E1	P147
											E1a	M132
												E1a1	M44
													E1a1a	A8384
														E1a1a1	Z17699
															E1a1a1a	Z17697
																E1a1a1a1	Z42331
																	E1a1a1a1a	Z17696
															E1a1a1b	A7710
														E1a1a2	Z17467
															E1a1a2a	Z31503
															E1a1a2b	SK569
												E1a2	CTS10713
													E1a2a	CTS246
														E1a2a1	CTS10935
															E1a2a1a	CTS3380
																E1a2a1a1	P110
																	E1a2a1a1a	CTS1792
																E1a2a1a2	Z15407
															E1a2a1b	L133.1
																E1a2a1b1	CTS736
																	E1a2a1b1a	CTS3764.1
														E1a2a2	Z5991
													E1a2b	Z5985.1
														E1a2b1	L94
															E1a2b1a	Z15064
																E1a2b1a1	Z5988
																	E1a2b1a1a	Z14979
																		E1a2b1a1a1	Z14889
																E1a2b1a2	Z5987
											E1b	P177
												E1b1	P2
													E1b1a	TSC0077541
														E1b1a1	M291
															E1b1a1a	P46
																E1b1a1a1	CTS1847
																	E1b1a1a1a	M4732
																		E1b1a1a1a1	CTS5078
																			E1b1a1a1a1a	M58
																				E1b1a1a1a1a1	Z6002
																				E1b1a1a1a1a2	CTS171
																			E1b1a1a1a1c	L485
																				E1b1a1a1a1c1	L514
																					E1b1a1a1a1c1a	CTS3046
																						E1b1a1a1a1c1a1	M191
																							E1b1a1a1a1c1a1a	P252
																								E1b1a1a1a1c1a1a3	CTS923
																									E1b1a1a1a1c1a1a3a	M4669
																										E1b1a1a1a1c1a1a3a1	CTS8030
																											E1b1a1a1a1c1a1a3a1a	M4437
																												E1b1a1a1a1c1a1a3a1a1	CTS9106
																													E1b1a1a1a1c1a1a3a1a1a	P116
																														E1b1a1a1a1c1a1a3a1a1a2	Z37839
																											E1b1a1a1a1c1a1a3a1b	Z5978
																												E1b1a1a1a1c1a1a3a1b1	Z5981
																												E1b1a1a1a1c1a1a3a1b2	CTS3642
																													E1b1a1a1a1c1a1a3a1b2a	Z5979
																														E1b1a1a1a1c1a1a3a1b2a1	Z5980
																											E1b1a1a1a1c1a1a3a1c	Z5982
																												E1b1a1a1a1c1a1a3a1c1	CTS67
																													E1b1a1a1a1c1a1a3a1c1a	CTS1871
																														E1b1a1a1a1c1a1a3a1c1a1	CTS2323
																													E1b1a1a1a1c1a1a3a1c1b	CTS11714
																														E1b1a1a1a1c1a1a3a1c1b1	Z1700
																														E1b1a1a1a1c1a1a3a1c1b2	Z22454
																															E1b1a1a1a1c1a1a3a1c1b2a	Z22570
																											E1b1a1a1a1c1a1a3a1d	CTS9806
																												E1b1a1a1a1c1a1a3a1d1	CTS8902
																													E1b1a1a1a1c1a1a3a1d1a	CTS1151
																													E1b1a1a1a1c1a1a3a1d1b	CTS9756
																														E1b1a1a1a1c1a1a3a1d1b1	Z1657
																															E1b1a1a1a1c1a1a3a1d1b1a	Z1655
																															E1b1a1a1a1c1a1a3a1d1b1b	CTS598
																																E1b1a1a1a1c1a1a3a1d1b1b1	Z1657
																																	E1b1a1a1a1c1a1a3a1d1b1b1a	Z1655
																														E1b1a1a1a1c1a1a3a1d1b2	CTS1338
																													E1b1a1a1a1c1a1a3a1d1c	CTS252
																													E1b1a1a1a1c1a1a3a1d1d	CTS89
																									E1b1a1a1a1c1a1a3b	CTS5299
																									E1b1a1a1a1c1a1a3c	CTS2754
																										E1b1a1a1a1c1a1a3c1	CTS5788
																											E1b1a1a1a1c1a1a3c1a	CTS6454
																											E1b1a1a1a1c1a1a3c1b	CTS1313
																												E1b1a1a1a1c1a1a3c1b1	Z5963
																													E1b1a1a1a1c1a1a3c1b1a	Z5964
																												E1b1a1a1a1c1a1a3c1b2	Z5965
																										E1b1a1a1a1c1a1a3c2	Z5966
																											E1b1a1a1a1c1a1a3c2a	Z5967
																											E1b1a1a1a1c1a1a3c2b	Z5968
																												E1b1a1a1a1c1a1a3c2b1	Z5969
																													E1b1a1a1a1c1a1a3c2b1a	Z5970
																														E1b1a1a1a1c1a1a3c2b1a1	Z22400
																															E1b1a1a1a1c1a1a3c2b1a1a	Z5972
																									E1b1a1a1a1c1a1a3d	CTS668
																										E1b1a1a1a1c1a1a3d4	CTS6508
																										E1b1a1a1a1c1a1a3d5	Z5974
																										E1b1a1a1a1c1a1a3d6	Z5975
																											E1b1a1a1a1c1a1a3d6a	V1419
																									E1b1a1a1a1c1a1a3e	Z22359
																									E1b1a1a1a1c1a1a3l	CTS553
																									E1b1a1a1a1c1a1a3m	Z22636
																								E1b1a1a1a1c1a1a4	Z5961
																							E1b1a1a1a1c1a1b	M4311
																								E1b1a1a1a1c1a1b1	CTS4662
																					E1b1a1a1a1c1b	L515
																						E1b1a1a1a1c1b1	CTS3393
																						E1b1a1a1a1c1b2	Z22239
																							E1b1a1a1a1c1b2a	M4273
																								E1b1a1a1a1c1b2a1	M4274
																									E1b1a1a1a1c1b2a1a	Z6027
																								E1b1a1a1a1c1b2a2	Z6025
																								E1b1a1a1a1c1b2a3	Z6023
																									E1b1a1a1a1c1b2a3a	Z6024
																				E1b1a1a1a1c2	CTS9883
																					E1b1a1a1a1c2a	Z6003
																						E1b1a1a1a1c2a1	Z6004
																					E1b1a1a1a1c2b	Z6005
																						E1b1a1a1a1c2b1	Z6006
																						E1b1a1a1a1c2b2	Z6007
																							E1b1a1a1a1c2b2a	Z16052
																					E1b1a1a1a1c2c	CTS3274
																						E1b1a1a1a1c2c1	CTS12004
																							E1b1a1a1a1c2c1a	CTS3725
																							E1b1a1a1a1c2c1b	Z6009
																						E1b1a1a1a1c2c2	Z6010
																							E1b1a1a1a1c2c2a	F2481
																						E1b1a1a1a1c2c3	Z6011
																							E1b1a1a1a1c2c3a	Z6012
																								E1b1a1a1a1c2c3a1	Z6013
																									E1b1a1a1a1c2c3a1a	Z6014
																									E1b1a1a1a1c2c3a1b	Z6015
																										E1b1a1a1a1c2c3a1b1	Z6016
																								E1b1a1a1a1c2c3a2	Z6017
																									E1b1a1a1a1c2c3a2a	Z6018
																								E1b1a1a1a1c2c3a3	Z6019
																							E1b1a1a1a1c2c3b	Z6020
																		E1b1a1a1a2	M4231
																			E1b1a1a1a2a	U175
																				E1b1a1a1a2a1	M4254
																					E1b1a1a1a2a1a	M4233
																						E1b1a1a1a2a1a1	Z21660
																						E1b1a1a1a2a1a2	CTS1433
																							E1b1a1a1a2a1a2a	CTS618
																						E1b1a1a1a2a1a3	M4257
																							E1b1a1a1a2a1a3a	CTS1974
																								E1b1a1a1a2a1a3a1	CTS9954
																									E1b1a1a1a2a1a3a1a	Z5947
																										E1b1a1a1a2a1a3a1a1	Z5948
																											E1b1a1a1a2a1a3a1a1a	Z5949
																								E1b1a1a1a2a1a3a2	CTS2198
																									E1b1a1a1a2a1a3a2a	F1823
																										E1b1a1a1a2a1a3a2a1	Z5942
																											E1b1a1a1a2a1a3a2a1a	Z5943
																												E1b1a1a1a2a1a3a2a1a1	Z8036
																												E1b1a1a1a2a1a3a2a1a2	Z5944
																												E1b1a1a1a2a1a3a2a1a3	Z5946
																												E1b1a1a1a2a1a3a2a1a4	Z42676
																							E1b1a1a1a2a1a3b	CTS6779
																								E1b1a1a1a2a1a3b1	CTS11
																									E1b1a1a1a2a1a3b1a	U290
																										E1b1a1a1a2a1a3b1a1	U181
																										E1b1a1a1a2a1a3b1a2	CTS99
																											E1b1a1a1a2a1a3b1a2a	M3912
																												E1b1a1a1a2a1a3b1a2a1	CTS127
																												E1b1a1a1a2a1a3b1a2a2	CTS1387
																													E1b1a1a1a2a1a3b1a2a2a	CTS907
																														E1b1a1a1a2a1a3b1a2a2a1	Z1761
																										E1b1a1a1a2a1a3b1a4	Z5950
																										E1b1a1a1a2a1a3b1a8	Z5954
																											E1b1a1a1a2a1a3b1a8a	CTS2184
																										E1b1a1a1a2a1a3b1a9	Z8038
																											E1b1a1a1a2a1a3b1a9a	Z5951
																										E1b1a1a1a2a1a3b1a10	CTS421
																										E1b1a1a1a2a1a3b1a11	Z5952
																										E1b1a1a1a2a1a3b1a12	Z5955
																											E1b1a1a1a2a1a3b1a12a	Z8037
																										E1b1a1a1a2a1a3b1a19	L651
																									E1b1a1a1a2a1a3b1d	M3905
																											E1b1a1a1a2a1a3b1d1a	Z5959
																											E1b1a1a1a2a1a3b1d1b	Z21511
																											E1b1a1a1a2a1a3b1d1c	CTS11328
																												E1b1a1a1a2a1a3b1d1c1	CTS3243
																													E1b1a1a1a2a1a3b1d1c1a	CTS10652
																														E1b1a1a1a2a1a3b1d1c1a1	CTS283
																															E1b1a1a1a2a1a3b1d1c1a1a	CTS1458
																								E1b1a1a1a2a1a3b3	CTS5768
																				E1b1a1a1a2a2	CTS1036
																					E1b1a1a1a2a2a	CTS6143
																						E1b1a1a1a2a2a1	Z21623
																						E1b1a1a1a2a2a2	CTS10239
																						E1b1a1a1a2a2a3	Z5940
																					E1b1a1a1a2a2b	M3940
															E1b1a1b	Z5994
																E1b1a1b1	Z22132
																	E1b1a1b1a	Z22149
																		E1b1a1b1a1	Z6000
																	E1b1a1b1b	Z5996
																	E1b1a1b1c	Z22205
																E1b1a1b2	Z36368
														E1b1a2	M329
													E1b1b	M215
														E1b1b1	M35.1
															E1b1b1a	L539
																E1b1b1a1	M78
																	E1b1b1a1a	CTS10890
																		E1b1b1a1a1	V12
																			E1b1b1a1a1b	V32
																				E1b1b1a1a1b1	CTS3282
																						E1b1b1a1a1b1a3	CTS1210
																			E1b1b1a1a1c	CTS693
																				E1b1b1a1a1c1	CTS3346
																					E1b1b1a1a1c1a	CTS1239
																						E1b1b1a1a1c1a1	Z5006
																				E1b1b1a1a1c2	CTS10132
																		E1b1b1a1a2	PF2272
																			E1b1b1a1a2a	PF2187
																				E1b1b1a1a2a1	CTS194
																	E1b1b1a1b	CTS4231
																		E1b1b1a1b1	L618
																			E1b1b1a1b1a	L142.1
																				E1b1b1a1b1a1	L17
																				E1b1b1a1b1a9	Z17293
																				E1b1b1a1b1a10	S2979
																					E1b1b1a1b1a10a	Z16659
																						E1b1b1a1b1a10a1	Y3183
																							E1b1b1a1b1a10a1a	Z16660
																								E1b1b1a1b1a10a1a2	Z36854
																						E1b1b1a1b1a10a2	Z21357
																					E1b1b1a1b1a10b	B409
																						E1b1b1a1b1a10b2	FGC11450
																		E1b1b1a1b2	CTS9547
																			E1b1b1a1b2a	CTS11457
																				E1b1b1a1b2a3	Z22640
															E1b1b1b	CTS1243
																E1b1b1b1	M310.1
																	E1b1b1b1a	M81
																		E1b1b1b1a1	M183
																			E1b1b1b1a1b	M5043
																			E1b1b1b1a1c	Y8829
																				E1b1b1b1a1c3	Z5010
																	E1b1b1b1b	PF2431
																		E1b1b1b1b1	FGC18905
																			E1b1b1b1b1b	PF2421
																E1b1b1b2	PF1961
																	E1b1b1b2a	PF1962
																		E1b1b1b2a1	M123
																				E1b1b1b2a1a1	CTS4483
																											E1b1b1b2a1a1a1a1b2a	Z36149
																											E1b1b1b2a1a1a1a1b2b	Y16781
																												E1b1b1b2a1a1a1a1b2b1	Y17226
																													E1b1b1b2a1a1a1a1b2b1a	Y17227
																														E1b1b1b2a1a1a1a1b2b1a1	Y18621
																															E1b1b1b2a1a1a1a1b2b1a1a	A16692
																															E1b1b1b2a1a1a1a1b2b1a1b	A12255
																															E1b1b1b2a1a1a1a1b2b1a1c	A20218
																															E1b1b1b2a1a1a1a1b2b1a1d	A15856
																														E1b1b1b2a1a1a1a1b2b1a2	Y17357
																												E1b1b1b2a1a1a1a1b2b2	Y145365
																									E1b1b1b2a1a1a1a1c	M136
																											E1b1b1b2a1a1a1a1f1b	PF6747
																												E1b1b1b2a1a1a1a1f1b1	BY1098
																													E1b1b1b2a1a1a1a1f1b1a	Z21018
																														E1b1b1b2a1a1a1a1f1b1a1	Z5001
																														E1b1b1b2a1a1a1a1f1b1a2	FGC15536
																				E1b1b1b2a1a2	M290
																				E1b1b1b2a1a4	L791
																	E1b1b1b2b	CTS10880
																		E1b1b1b2b2	CTS1177
																				E1b1b1b2b2a1	M293
																							E1b1b1b2b2a1a1a	Z865
															E1b1b1c	V6
															E1b1b1d	V92
														E1b1b2	M281
												E1b2	P75
										E2	M75
											E2b	M98
												E2b1	M54
													E2b1a	M85
														E2b1a1	M200
															E2b1a1c	CTS2629
															E2b1a1d	CTS388
												E2b2	CTS1441
								CF	P143
									C	M130
										C1	F3393
											C1a	CTS11043
												C1a1	M8
													C1a1a	P121
														C1a1a1	CTS9336
															C1a1a1a	Z7180
																C1a1a1a1	CTS6678
															C1a1a1b	Z1356
																C1a1a1b2	Z40544
												C1a2	V20
													C1a2a	V182
														C1a2a1	V222
															C1a2a1a	Y12152
																C1a2a1a1	BY1117
																	C1a2a1a1a	Z31799
																C1a2a1a2	Y12157
																	C1a2a1a2a	BY44466
															C1a2a1b	BY67541
																C1a2a1b1	BY98578
														C1a2a2	Y11325
													C1a2b	Z38888
														C1a2b1	BY22666
															C1a2b1a	Z44576
																C1a2b1a1	Y140262
															C1a2b1b	Z44526
																C1a2b1b1	PH188.2
																C1a2b1b2	F18203
											C1b	F1370
												C1b1	AM00694
													C1b1a	B66
														C1b1a1	SK991
																C1b1a1a1	P92
																	C1b1a1a1a	Z5895
																		C1b1a1a1a1	K159
																			C1b1a1a1a1a	K96
																				C1b1a1a1a1a1	K42
																					C1b1a1a1a1a1a	Z5896
																						C1b1a1a1a1a1a1	K193
																							C1b1a1a1a1a1a1a	Z12522
																								C1b1a1a1a1a1a1a1	K466
																									C1b1a1a1a1a1a1a1a	K470
																										C1b1a1a1a1a1a1a1a1	Z5898
																C1b1a1a2	Z5899
																	C1b1a1a2a	Z12338
																		C1b1a1a2a2	SK994
															C1b1a1b	Z16582
														C1b1a2	AM00847
															C1b1a2a	B67
															C1b1a2b	AM00594
																C1b1a2b1	AM00627
																	C1b1a2b1a	B465
												C1b2	B477
													C1b2a	M38
														C1b2a1	M208
																C1b2a1b1	P54
															C1b2a1c	F21693
																C1b2a1c1	B460
																	C1b2a1c1b	B461
																C1b2a1c2	FT12310
																	C1b2a1c2a	B459
																	C1b2a1c2b	FT12484
																		C1b2a1c2b1	FT13521
																			C1b2a1c2b1a	FT13115
																		C1b2a1c2b2	FT134718
																				C1b2a1c2b2a1	FT133627
																		C1b2a1c2b3	FT232409
																			C1b2a1c2b3a	FT245596
														C1b2a2	FT71404
															C1b2a2b	F2642.3
													C1b2b	M347
														C1b2b1	M210
										C2	M217
											C2a	L1373
												C2a1	F3447
													C2a1a	F1699
														C2a1a1	F3918
															C2a1a1a	P39
																C2a1a1a1	BY1360
																	C2a1a1a1a	Z30538
																		C2a1a1a1a1	Z30750
																			C2a1a1a1a1a	Z30754
																		C2a1a1a1a2	BY736
																			C2a1a1a1a2a	Z30765
																			C2a1a1a1a2c	BY22870
																	C2a1a1a1b	Z38940
																C2a1a1a2	Z38874
															C2a1a1b	FGC28881.2
																C2a1a1b1	F1756
																	C2a1a1b1a	F3830
																		C2a1a1b1a1	Z32845
																		C2a1a1b1a2	F10085
																			C2a1a1b1a2a	F9721
																				C2a1a1b1a2a1	FGC28846
																					C2a1a1b1a2a1b	FGC31362
																				C2a1a1b1a2a2	BY186297
																	C2a1a1b1b	Y10420
																		C2a1a1b1b1	Y11606
														C2a1a2	M48
															C2a1a2a	M86
																	C2a1a2a1a	B469
																		C2a1a2a1a1	B87
																			C2a1a2a1a1a	B89
																	C2a1a2a1b	B80
																		C2a1a2a1b1	B471
																		C2a1a2a1b2	B81
																			C2a1a2a1b2b	B82
																	C2a1a2a2a	F5485
																		C2a1a2a2a1	Y15844
																			C2a1a2a2a1a	Y15552
																				C2a1a2a2a1a1	BY18743
																				C2a1a2a2a1a2	MPB421
																				C2a1a2a2a1a3	Y17342
																		C2a1a2a2a2	SK1066
																		C2a1a2a2a3	F12970
																			C2a1a2a2a3a	F9936
																				C2a1a2a2a3a1	F9766
																	C2a1a2a2b	F8472
															C2a1a2b	B90
																C2a1a2b1	B91
																	C2a1a2b1b	B92
														C2a1a3	M504
															C2a1a3a	F3796
																C2a1a3a1	F5481
																	C2a1a3a1c	FGC16217
																	C2a1a3a1e	Y12782
																		C2a1a3a1e2	F5483
																		C2a1a3a1e3	BY18686
																	C2a1a3a1f	F10091
																	C2a1a3a1i	BY182928
																C2a1a3a6	SK1072
																	C2a1a3a6a	Y148768
																		C2a1a3a6a1	FT37001
																	C2a1a3a6b	F12663
														C2a1a4	F9992
															C2a1a4a	Z22425
																	C2a1a4a1a	BY56585
													C2a1b	Z31698
												C2a2	BY63635
													C2a2b	BY67935
											C2b	F1067
												C2b1	F2613
													C2b1a	Z1300
														C2b1a1	CTS2657
															C2b1a1a	CTS11990
																C2b1a1a1	CTS8579
																	C2b1a1a1a	M407
																		C2b1a1a1a1	F3850
																			C2b1a1a1a1a	B97
																				C2b1a1a1a1a5	B98
																			C2b1a1a1a1c	B99
																		C2b1a1a1a2	Z45401
																C2b1a1a2	Z31664
																	C2b1a1a2a	Z31665
															C2b1a1b	A14895. A14896
																C2b1a1b1	A14909
														C2b1a2	K700
															C2b1a2a	F1319
																C2b1a2a1	F3777
																	C2b1a2a1a	F9966. F11237
																		C2b1a2a1a1	F3748
																		C2b1a2a1a2	BY53999
																			C2b1a2a1a2a	FT40406
																C2b1a2a2	F10056
																	C2b1a2a2b	F3555
																		C2b1a2a2b1	F10036
																			C2b1a2a2b1a	F13864
																				C2b1a2a2b1a2	F13136
																					C2b1a2a2b1a2b	BY96205
																						C2b1a2a2b1a2b1	MF4473
															C2b1a2b	CTS3385
																C2b1a2b2	FGC45548
													C2b1b	F845
														C2b1b1	CTS10923
															C2b1b1a	K511
																C2b1b1a1	K516
															C2b1b1b	F5477
																C2b1b1b2	SK1038
																	C2b1b1b2a	MF1017
																		C2b1b1b2a1	MF1020
															C2b1b1d	FGC39587
																C2b1b1d2	Z43981
															C2b1b1g	Z45293
														C2b1b2	MF2091
												C2b2	CTS4660
									F	M89
										F	M89
											F1	P91
											F2	M427
											F3	M481
											GHIJK	F1329
										GHIJK	F1329
											G	M201
												G1	M342
													G1a	CTS11562
														G1a1	GG349
															G1a1a	BY1124
																G1a1a1	Z3353
																	G1a1a1a	L1324
																		G1a1a1a1	Y12260
																			G1a1a1a1a	L201
																			G1a1a1a1d	Z44585
																			G1a1a1a1e	FT107951
																				G1a1a1a1e2	FT107429
																		G1a1a1a2	L1323
																	G1a1a1b	GG362
																		G1a1a1b1	Z26332
																			G1a1a1b1a	FT101799
																		G1a1a1b2	FT17156
																			G1a1a1b2b	FT19393
																				G1a1a1b2b2	FT19860
																G1a1a2	GG313
																	G1a1a2a	GG167
																		G1a1a2a1	GG169
																			G1a1a2a1a	GG164
																	G1a1a2b	GG272
																		G1a1a2b1	GG270
																			G1a1a2b1a	GG265
																				G1a1a2b1a2	Y134359
																			G1a1a2b1c	BY175028
																				G1a1a2b1c2	BY175888
																						G1a1a2b1c2b1	BY176232
																		G1a1a2b2	BY94396
																			G1a1a2b2a	BY60348
																				G1a1a2b2a1	BY157285
																					G1a1a2b2a1a	BY73830
																				G1a1a2b2a2	Z45949
														G1a2	F2885
															G1a2a	BY21302
																G1a2a1	BY95017
																G1a2a2	Z45001
																	G1a2a2a	BY44900
																		G1a2a2a1	Z45002
													G1b	L830
														G1b1	Z17874
															G1b1a	Z18606
																G1b1a1	BY56952
															G1b1b	Z30744
												G2	P287
													G2a	P15
														G2a1	FGC7535
														[L293  used in several studies to represent the FGC7535 subgroup was withdrawn due to problems]	[L293  used in several studies to represent the FGC7535 subgroup was withdrawn due to problems]
															G2a1a	FGC595
																G2a1a1	FGC776
																	G2a1a1a	FGC693
																		G2a1a1a1	FGC715
																			G2a1a1a1a	Z6635
																				G2a1a1a1a1	Z6638
																					G2a1a1a1a1a	FGC694
																						G2a1a1a1a1a1	FGC701
																							G2a1a1a1a1a1a	FGC766
																								G2a1a1a1a1a1a1	Z7958
																									G2a1a1a1a1a1a1a	FGC750
																										G2a1a1a1a1a1a1a1	FGC724
																											G2a1a1a1a1a1a1a1a	Z7961
																												G2a1a1a1a1a1a1a1a1	Z7944
																													G2a1a1a1a1a1a1a1a1a	FGC3762.2
																														G2a1a1a1a1a1a1a1a1a2	FGC3750.2
																															G2a1a1a1a1a1a1a1a1a2a	FGC3787
																																G2a1a1a1a1a1a1a1a1a2a1	FGC7909
																																	G2a1a1a1a1a1a1a1a1a2a1b	FT209164
																																		G2a1a1a1a1a1a1a1a1a2a1b1	FT286824
																														G2a1a1a1a1a1a1a1a1a3	FT61413
																												G2a1a1a1a1a1a1a1a2	FGC719
																													G2a1a1a1a1a1a1a1a2a	FGC668
																														G2a1a1a1a1a1a1a1a2a1	A199
																														G2a1a1a1a1a1a1a1a2a3	FGC690
																													G2a1a1a1a1a1a1a1a2b	FT185956
																														G2a1a1a1a1a1a1a1a2b1	FT183164
																												G2a1a1a1a1a1a1a1a3	Z40550
																											G2a1a1a1a1a1a1a1b	FT70700
																										G2a1a1a1a1a1a1a2	Z31459
																											G2a1a1a1a1a1a1a2a	Y36736
																												G2a1a1a1a1a1a1a2a1	Z31461
																													G2a1a1a1a1a1a1a2a1a	BY21462
																														G2a1a1a1a1a1a1a2a1a1	Y167948
																															G2a1a1a1a1a1a1a2a1a1a	Z31460
																														G2a1a1a1a1a1a1a2a1a2	Z45052
																															G2a1a1a1a1a1a1a2a1a2a	Z45048
																												G2a1a1a1a1a1a1a2a2	Y97597
																											G2a1a1a1a1a1a1a2c	BY164577
																								G2a1a1a1a1a1a2	FGC1047
																									G2a1a1a1a1a1a2a	FGC1048
																					G2a1a1a1a1b	FT19842
																			G2a1a1a1b	FGC1160
																				G2a1a1a1b1	Z6673
																					G2a1a1a1b1a	FGC1107
																				G2a1a1a1b3	BY21464
																					G2a1a1a1b3a	Z29313
																						G2a1a1a1b3a1	Z29311
																		G2a1a1a2	Z31455
																			G2a1a1a2b	Z25119.2
																		G2a1a1a3	FGC65315
																			G2a1a1a3a	BY201496
																G2a1a2	Z17774
																	G2a1a2a	Z17775
																		G2a1a2a1	Z31475
																			G2a1a2a1a	Z39928
																G2a1a3	Z31464
																	G2a1a3a	Z35447
																	G2a1a3b	FGC50490
															G2a1b	BY65109
																G2a1b2	BY198464
																	G2a1b2a	FT66003
														G2a2	CTS4367
															G2a2a	PF3147
																G2a2a1	PF3148
																	G2a2a1a	PF3177
																		G2a2a1a1	FGC6669
																			G2a2a1a1a	FGC58131
																				G2a2a1a1a1	M286
																					G2a2a1a1a1a	FGC58189
																						G2a2a1a1a1a2	FGC59473
																						G2a2a1a1a1b1	Z43083
																				G2a2a1a1b1	FGC6618
																					G2a2a1a1b1a	FGC6629
																		G2a2a1a2	L91
																			G2a2a1a2a	FGC7739
																				G2a2a1a2a1	PF3239
																					G2a2a1a2a1a	L166
																						G2a2a1a2a1a1	FGC5672
																							G2a2a1a2a1a1a	Z6211
																								G2a2a1a2a1a1a1	FGC5697
																								G2a2a1a2a1a1b1	Z31442
																					G2a2a1a2a1b	Z6802
																								G2a2a1a2a1b1a2	Z31447 BY180234
																				G2a2a1a2a2	FGC2315
																					G2a2a1a2a2a	Z6773
																						G2a2a1a2a2a1	Y140837
																							G2a2a1a2a2a1a	Z31429
																								G2a2a1a2a2a1a1	Z31428
																									G2a2a1a2a2a1a1a	BY81231
																								G2a2a1a2a2a1a2	Y140821
																							G2a2a1a2a2a1b	FT99363
																					G2a2a1a2a2b	FGC2308
																			G2a2a1a2b	Z42562
																				G2a2a1a2b1	Z42570
																					G2a2a1a2b1b	Z42554
																				G2a2a1a2b2	BY45264
																		G2a2a1a3	FGC34625
																		G2a2a1a5	BY182515
																	G2a2a1b	FGC34451
																		G2a2a1b2	FGC35165
																G2a2a2	Z36520
																	G2a2a2a	Z45970
																		G2a2a2a1	Z31724
																	G2a2a2b	Z36525
																		G2a2a2b1	BY37102
																			G2a2a2b1a	Y90753
																				G2a2a2b1a1	Y83453
																				G2a2a2b1a3	BY193093
																		G2a2a2b3	Z45700
																			G2a2a2b3a	Z45690
																							G2a2a2b3a1a1c	FT19778
															G2a2b	L30
																G2a2b1	M406
																	G2a2b1a	FGC5089
																		G2a2b1a1	FGC5081
																			G2a2b1a1a	FGC5091
																				G2a2b1a1a1	L14
																					G2a2b1a1a1a	FGC5185
																						G2a2b1a1a1a1	Z17085
																							G2a2b1a1a1a1a	Z17083
																						G2a2b1a1a1a2	Y92117
																							G2a2b1a1a1a2b	BY91399
																								G2a2b1a1a1a2b2	FT44209
																									G2a2b1a1a1a2b2a	FT40739
																						G2a2b1a1a1a7	Y38327
																						G2a2b1a1a1a9	FT45116
																							G2a2b1a1a1a9a	FT167155
																					G2a2b1a1a1b	Z45043
																						G2a2b1a1a1b2	Y128488
																							G2a2b1a1a1b2a	Z45580
																								G2a2b1a1a1b2a1	Z45581
																						G2a2b1a1a1b3	FT115489
																			G2a2b1a1b	Z17887
																				G2a2b1a1b1	Z17886
																					G2a2b1a1b1a	Z37368
																						G2a2b1a1b1a2	Z37375
																						G2a2b1a1b1a3	Z39512
																				G2a2b1a1b2	PH942
																					G2a2b1a1b2a	L645
																						G2a2b1a1b2a1	Z42373
																								G2a2b1a1b2a1a1	Z21254.2
																									G2a2b1a1b2a1a1a	Z42385
																									G2a2b1a1b2a1a1b	BY93324
																						G2a2b1a1b2a2	PH1471
																						G2a2b1a1b2a3	BY61547
																							G2a2b1a1b2a3a	FT108136
																			G2a2b1a1e	S9350
																				G2a2b1a1e1	Y24772
																					G2a2b1a1e1a	S10406
																					G2a2b1a1e1b	Y24774
																						G2a2b1a1e1b1	Y24444
																				G2a2b1a1e2	BY172388
																		G2a2b1a2	M3302
																			G2a2b1a2a	M3422
																				G2a2b1a2a1	M3240
																			G2a2b1a2b	Z31412
																				G2a2b1a2b1	BY21758
																					G2a2b1a2b1a	BY21760
																						G2a2b1a2b1a1	M5843
																						G2a2b1a2b1a2	BY69342
																							G2a2b1a2b1a2a	BY53866
																							G2a2b1a2b1a2b	BY177020.2
																			G2a2b1a2c	FT56736
																	G2a2b1b	PF3293
																		G2a2b1b1	S25020
																			G2a2b1b1a	PF3296
																				G2a2b1b1a1	Z30795
																					G2a2b1b1a1a	Z30797
																						G2a2b1b1a1a1	Z30793
																							G2a2b1b1a1a1c	FT75972
																								G2a2b1b1a1a1c1	BY68816
																					G2a2b1b1a1b	S11415
																						G2a2b1b1a1b1	Y82047
																							G2a2b1b1a1b1a	Z45959
																							G2a2b1b1a1b1c	BY75675
																								G2a2b1b1a1b1c2	FT29873
																							G2a2b1b1a1b1d	BY187298
																								G2a2b1b1a1b1d1	BY186341
																									G2a2b1b1a1b1d1a	BY197732
																										G2a2b1b1a1b1d1a2	BY197830
																						G2a2b1b1a1b3	BY104772
																							G2a2b1b1a1b3a	BY83406
																						G2a2b1b1a1b4	FT144617
																							G2a2b1b1a1b4a	FT259551
																						G2a2b1b1a1b5	BY169001
																					G2a2b1b1a1c	Z6354
																						G2a2b1b1a1c1	Z6352
																				G2a2b1b1a2	Y32612
																					G2a2b1b1a2a	Y32613
																						G2a2b1b1a2a1	Y94775
																					G2a2b1b1a2b	Z28113.2
																	G2a2b1d	CTS8450
																G2a2b2	CTS2488
																	G2a2b2a	P303
																		G2a2b2a1	L140
																			G2a2b2a1a	PF3346
																				G2a2b2a1a1	PF3345
																					G2a2b2a1a1a	FGC7568
																						G2a2b2a1a1a1	U1
																							G2a2b2a1a1a1a	L13
																								G2a2b2a1a1a1a1	CTS9909
																									G2a2b2a1a1a1a1a	FGC965
																										G2a2b2a1a1a1a1a1	FGC998
																											G2a2b2a1a1a1a1a1a	FGC1016
																												G2a2b2a1a1a1a1a1a1	S9751
																													G2a2b2a1a1a1a1a1a1a	Z29424
																														G2a2b2a1a1a1a1a1a1a1	L1263
																															G2a2b2a1a1a1a1a1a1a1a	Z38846
																																G2a2b2a1a1a1a1a1a1a1a1	BY48807
																																	G2a2b2a1a1a1a1a1a1a1a1a	Y19318
																																		G2a2b2a1a1a1a1a1a1a1a1a1	Z38853
																																		G2a2b2a1a1a1a1a1a1a1a1a3	Z43288
																																		G2a2b2a1a1a1a1a1a1a1a1a4	Z43941
																																		G2a2b2a1a1a1a1a1a1a1a1a8	FT34806
																																		G2a2b2a1a1a1a1a1a1a1a1a10	BY37125
																														G2a2b2a1a1a1a1a1a1a2	Z30777
																															G2a2b2a1a1a1a1a1a1a2a	Z30781
																																G2a2b2a1a1a1a1a1a1a2a1	Z38455
																														G2a2b2a1a1a1a1a1a1a5	Z39053
																														G2a2b2a1a1a1a1a1a1a8	Z43944
																															G2a2b2a1a1a1a1a1a1a8b	BY46240
																														G2a2b2a1a1a1a1a1a1a9	Z44047
																														G2a2b2a1a1a1a1a1a1a11	Z45562
																															G2a2b2a1a1a1a1a1a1a11b	FT24243
																														G2a2b2a1a1a1a1a1a1a16	BY45246
																														G2a2b2a1a1a1a1a1a1a17	FT167239
																												G2a2b2a1a1a1a1a1a2	Z39373
																													G2a2b2a1a1a1a1a1a2a	Y45819
																														G2a2b2a1a1a1a1a1a2a1	Y57313
																											G2a2b2a1a1a1a1a1b	Z42471
																												G2a2b2a1a1a1a1a1b1	Z42474
																										G2a2b2a1a1a1a1a2	Z2017
																												G2a2b2a1a1a1a1a2a1	Z30718
																												G2a2b2a1a1a1a1a2a2	Z31377
																													G2a2b2a1a1a1a1a2a2a	FT213870
																									G2a2b2a1a1a1a1b	Z44912
																										G2a2b2a1a1a1a1b1	Z30831
																									G2a2b2a1a1a1a1c	S20738
																										G2a2b2a1a1a1a1c2	BY67743
																							G2a2b2a1a1a1b	L1266
																								G2a2b2a1a1a1b1	L1264
																									G2a2b2a1a1a1b1a	Z44222
																										G2a2b2a1a1a1b1a1	L654.2
																										G2a2b2a1a1a1b1a2	Y32922
																											G2a2b2a1a1a1b1a2a	Y32599
																												G2a2b2a1a1a1b1a2a1	Y32597
																												G2a2b2a1a1a1b1a2a2	FT8464
																											G2a2b2a1a1a1b1a2b	FT49803
																												G2a2b2a1a1a1b1a2b1	FT49834
																										G2a2b2a1a1a1b1a3	V7991.2
																											G2a2b2a1a1a1b1a3a	Y112447
																												G2a2b2a1a1a1b1a3a1	FT13327
																													G2a2b2a1a1a1b1a3a1a	FT12699
																													G2a2b2a1a1a1b1a3a1b	FT69519
																														G2a2b2a1a1a1b1a3a1b1	FT69517
																									G2a2b2a1a1a1b1b	FGC21495
																										G2a2b2a1a1a1b1b1	S9409
																											G2a2b2a1a1a1b1b1b	Y142068
																													G2a2b2a1a1a1b1b1b1a1	BY166419
																													G2a2b2a1a1a1b1b1b1a2	FT238907
																												G2a2b2a1a1a1b1b1b2	FT69558
																													G2a2b2a1a1a1b1b1b2a	BY54724
																										G2a2b2a1a1a1b1b2	Y30992
																											G2a2b2a1a1a1b1b2a	Z30715
																									G2a2b2a1a1a1b1c	BY109806
																										G2a2b2a1a1a1b1c1	BY79584
																										G2a2b2a1a1a1b1c2	BY4191
																										G2a2b2a1a1a1b1c3	FT145691
																										G2a2b2a1a1a1b1c4	FT178740
																								G2a2b2a1a1a1b2	PH1780
																									G2a2b2a1a1a1b2a	PH311.1
																										G2a2b2a1a1a1b2a4	Y128514
																							G2a2b2a1a1a1d	SK1154
																					G2a2b2a1a1b	L497
																						G2a2b2a1a1b1	CTS9737
																							G2a2b2a1a1b1a	CTS11352
																								G2a2b2a1a1b1a1	Z725
																									G2a2b2a1a1b1a1a	AMM042
																										G2a2b2a1a1b1a1a1	L43
																											G2a2b2a1a1b1a1a1a	L42
																												G2a2b2a1a1b1a1a1a1	Y11074
																													G2a2b2a1a1b1a1a1a1a	Y11076.2
																														G2a2b2a1a1b1a1a1a1a1	F1300.2
																															G2a2b2a1a1b1a1a1a1a1a	Z39505
																																G2a2b2a1a1b1a1a1a1a1a1	Z39504
																														G2a2b2a1a1b1a1a1a1a2	CTS7357.1
																															G2a2b2a1a1b1a1a1a1a2b	CTS3647.2
																																G2a2b2a1a1b1a1a1a1a2b2	Z31318
																														G2a2b2a1a1b1a1a1a1a3	Z31316
																															G2a2b2a1a1b1a1a1a1a3b	BY64222
																														G2a2b2a1a1b1a1a1a1a6	Z45511
																														G2a2b2a1a1b1a1a1a1a7	FGC57326
																													G2a2b2a1a1b1a1a1a1b	YSC0000033
																														G2a2b2a1a1b1a1a1a1b1	Z39501
																															G2a2b2a1a1b1a1a1a1b1b	Y128028
																																G2a2b2a1a1b1a1a1a1b1b2	BY195513
																														G2a2b2a1a1b1a1a1a1b2	M7109.2
																														G2a2b2a1a1b1a1a1a1b4	Z44675
																															G2a2b2a1a1b1a1a1a1b4a	Y98860
																																G2a2b2a1a1b1a1a1a1b4a1	Z44672
																																	G2a2b2a1a1b1a1a1a1b4a1b	Y88703
																																	G2a2b2a1a1b1a1a1a1b4a1c	Y138634
																															G2a2b2a1a1b1a1a1a1b4b	FT19044
																													G2a2b2a1a1b1a1a1a1c	Z45255
																														G2a2b2a1a1b1a1a1a1c2	Z41181
																														G2a2b2a1a1b1a1a1a1c3	Z45508
																															G2a2b2a1a1b1a1a1a1c3a	FGC34077
																															G2a2b2a1a1b1a1a1a1c3c	BY154149
																														G2a2b2a1a1b1a1a1a1c5	BY95079
																												G2a2b2a1a1b1a1a1a2	Z40854
																													G2a2b2a1a1b1a1a1a2a	Y65311.2
																														G2a2b2a1a1b1a1a1a2a2	BY140356
																															G2a2b2a1a1b1a1a1a2a2b	BY173771
																											G2a2b2a1a1b1a1a1b	Z31335
																												G2a2b2a1a1b1a1a1b1	Z31333
																										G2a2b2a1a1b1a1a2	CTS6796
																											G2a2b2a1a1b1a1a2a	CTS4803
																												G2a2b2a1a1b1a1a2a1	S2808
																													G2a2b2a1a1b1a1a2a1a	BY46757
																														G2a2b2a1a1b1a1a2a1a1	CTS10391
																															G2a2b2a1a1b1a1a2a1a1a	Z30787
																															G2a2b2a1a1b1a1a2a1a1b	S18765
																																G2a2b2a1a1b1a1a2a1a1b1	Z41143
																																	G2a2b2a1a1b1a1a2a1a1b1a	Z41144
																																		G2a2b2a1a1b1a1a2a1a1b1a2	BY55158
																																				G2a2b2a1a1b1a1a2a1a1b1a2b1	FT25136
																																					G2a2b2a1a1b1a1a2a1a1b1a2b1b	BY4514
																																						G2a2b2a1a1b1a1a2a1a1b1a2b1b1	FT19708
																																						G2a2b2a1a1b1a1a2a1a2a1b2b1b2	FT32514
																																	G2a2b2a1a1b1a1a2a1a1b1d	BY101101
																																		G2a2b2a1a1b1a1a2a1a1b1d1	Z45470
																																	G2a2b2a1a1b1a1a2a1a1b1e	BY87342
																																G2a2b2a1a1b1a1a2a1a1b2	Y106451
																																	G2a2b2a1a1b1a1a2a1a1b2a	BY87464
																																G2a2b2a1a1b1a1a2a1a1b3	BY114181
																															G2a2b2a1a1b1a1a2a1a1c	Z30790
																																G2a2b2a1a1b1a1a2a1a1c1	Z30791
																																G2a2b2a1a1b1a1a2a1a1c2	FGC34989
																														G2a2b2a1a1b1a1a2a1a2	FT167696
																													G2a2b2a1a1b1a1a2a1b	FGC8304
																														G2a2b2a1a1b1a1a2a1b1	FGC8303
																															G2a2b2a1a1b1a1a2a1b1a	S2795
																																G2a2b2a1a1b1a1a2a1b1a1	S2784
																																G2a2b2a1a1b1a1a2a1b1a2	Y103893
																															G2a2b2a1a1b1a1a2a1b1b	FGC8314
																																G2a2b2a1a1b1a1a2a1b1b1	FGC8354
																																	G2a2b2a1a1b1a1a2a1b1b1a	Z40548
																																	G2a2b2a1a1b1a1a2a1b1b1b	Z45516
																																G2a2b2a1a1b1a1a2a1b1b2	Z40765
																																	G2a2b2a1a1b1a1a2a1b1b2a	Z44626
																																	G2a2b2a1a1b1a1a2a1b1b2b	FT64136
																																G2a2b2a1a1b1a1a2a1b1b3	Z45548
																														G2a2b2a1a1b1a1a2a1b2	FGC34826
																															G2a2b2a1a1b1a1a2a1b2a	FGC34823
																																G2a2b2a1a1b1a1a2a1b2a1	FGC34846
																																	G2a2b2a1a1b1a1a2a1b2a1a	Y13110
																																		G2a2b2a1a1b1a1a2a1b2a1a1	Z30729
																																			G2a2b2a1a1b1a1a2a1b2a1a1a	Y16788
																																				G2a2b2a1a1b1a1a2a1b2a1a1a3	Y49262
																																			G2a2b2a1a1b1a1a2a1b2a1a1b	BY74723
																																				G2a2b2a1a1b1a1a2a1b2a1a1b2	BY28078
																																			G2a2b2a1a1b1a1a2a1b2a1a1c	Y30009
																																	G2a2b2a1a1b1a1a2a1b2a1b	FGC34822
																																		G2a2b2a1a1b1a1a2a1b2a1b1	FGC34830
																																	G2a2b2a1a1b1a1a2a1b2a1c	BY15210
																																G2a2b2a1a1b1a1a2a1b2a2	SK1149
																																G2a2b2a1a1b1a1a2a1b2a3	BY28015
																																	G2a2b2a1a1b1a1a2a1b2a3a	BY100747
																																		G2a2b2a1a1b1a1a2a1b2a3a2	PH4671
																																G2a2b2a1a1b1a1a2a1b2a4	Z17780
																																	G2a2b2a1a1b1a1a2a1b2a4a	Z45537
																																G2a2b2a1a1b1a1a2a1b2a5	Z46349
																															G2a2b2a1a1b1a1a2a1b2c	BY27991
																																G2a2b2a1a1b1a1a2a1b2c1	BY27998
																																G2a2b2a1a1b1a1a2a1b2c3	FT57386
																													G2a2b2a1a1b1a1a2a1c	FGC14522
																														G2a2b2a1a1b1a1a2a1c1	FGC14523
																												G2a2b2a1a1b1a1a2a2	Z6150
																													G2a2b2a1a1b1a1a2a2a	Z30771
																														G2a2b2a1a1b1a1a2a2a1	F3484.2
																															G2a2b2a1a1b1a1a2a2a1a	Z17787
																														G2a2b2a1a1b1a1a2a2a2	BY28016
																													G2a2b2a1a1b1a1a2a2b	FT100320
																														G2a2b2a1a1b1a1a2a2b1	Z44661
																															G2a2b2a1a1b1a1a2a2b1a	Z44658
																																G2a2b2a1a1b1a1a2a2b1a1	Y93363
																														G2a2b2a1a1b1a1a2a2b2	BY117124
																															G2a2b2a1a1b1a1a2a2b2a	BY138628
																															G2a2b2a1a1b1a1a2a2b2b	FT100736
																													G2a2b2a1a1b1a1a2a2c	FT6731
																													G2a2b2a1a1b1a1a2a2d	Z40694
																												G2a2b2a1a1b1a1a2a3	BY102306
																												G2a2b2a1a1b1a1a2a4	Z42514
																													G2a2b2a1a1b1a1a2a4a	Z31329
																													G2a2b2a1a1b1a1a2a4c	BY37149
																												G2a2b2a1a1b1a1a2a5	FGC60011
																												G2a2b2a1a1b1a1a2a12	S15656
																													G2a2b2a1a1b1a1a2a12b	FGC36009
																												G2a2b2a1a1b1a1a2a14	BY199514
																													G2a2b2a1a1b1a1a2a14a	FT126341
																											G2a2b2a1a1b1a1a2b	Z16775
																												G2a2b2a1a1b1a1a2b1	Y14684
																													G2a2b2a1a1b1a1a2b1a	Z16770
																														G2a2b2a1a1b1a1a2b1a1	FGC21278
																															G2a2b2a1a1b1a1a2b1a1a	Y14683
																																G2a2b2a1a1b1a1a2b1a1a1	Z34079
																																	G2a2b2a1a1b1a1a2b1a1a1a	Z34083
																																		G2a2b2a1a1b1a1a2b1a1a1a1	Z34087
																																	G2a2b2a1a1b1a1a2b1a1a1b	Z31343
																																G2a2b2a1a1b1a1a2b1a1a2	Z31842
																																	G2a2b2a1a1b1a1a2b1a1a2a	Z37851
																																G2a2b2a1a1b1a1a2b1a1a3	Y68656
																																	G2a2b2a1a1b1a1a2b1a1a3a	Z31338
																																G2a2b2a1a1b1a1a2b1a1a6	BY62693
																																	G2a2b2a1a1b1a1a2b1a1a6a	Z42370
																																G2a2b2a1a1b1a1a2b1a1a8	FT177151
																														G2a2b2a1a1b1a1a2b1a2	BY179461
																													G2a2b2a1a1b1a1a2b1b	Z41650
																														G2a2b2a1a1b1a1a2b1b1	FT230034
																															G2a2b2a1a1b1a1a2b1b1a	Z41649
																														G2a2b2a1a1b1a1a2b1b2	BY123060
																														G2a2b2a1a1b1a1a2b1b4	FT4953
																												G2a2b2a1a1b1a1a2b2	Z27567
																													G2a2b2a1a1b1a1a2b2a	Z39670
																														G2a2b2a1a1b1a1a2b2a1	Z39674
																															G2a2b2a1a1b1a1a2b2a1a	Z39673
																																G2a2b2a1a1b1a1a2b2a1a1	Z39688
																																	G2a2b2a1a1b1a1a2b2a1a1a	BY45702
																																G2a2b2a1a1b1a1a2b2a1a2	Z39689
																																G2a2b2a1a1b1a1a2b2a1a3	Z39694
																																G2a2b2a1a1b1a1a2b2a1a5	Z45631
																														G2a2b2a1a1b1a1a2b2a2	Z40538
																											G2a2b2a1a1b1a1a2c	Z36217
																												G2a2b2a1a1b1a1a2c1	FT115823
																													G2a2b2a1a1b1a1a2c1a1	FT115715
																										G2a2b2a1a1b1a1a3	Z30724
																											G2a2b2a1a1b1a1a3a	Z30725
																												G2a2b2a1a1b1a1a3a1	BY28079
																										G2a2b2a1a1b1a1a4	Z44678. BY120911
																											G2a2b2a1a1b1a1a4a	Z44677
																									G2a2b2a1a1b1a1b	FGC477
																										G2a2b2a1a1b1a1b1	FGC489
																											G2a2b2a1a1b1a1b1a	FGC490
																												G2a2b2a1a1b1a1b1a1	FGC486
																												G2a2b2a1a1b1a1b1a2	Z40857
																													G2a2b2a1a1b1a1b1a2a	FT21948
																														G2a2b2a1a1b1a1b1a2a1	Z40858
																															G2a2b2a1a1b1a1b1a2a1a	FT68913
																														G2a2b2a1a1b1a1b1a2a2	Y132507
																															G2a2b2a1a1b1a1b1a2a2a	BY94332
																														G2a2b2a1a1b1a1b1a2a3	FT139285
																													G2a2b2a1a1b1a1b1a2b	Y132505
																														G2a2b2a1a1b1a1b1a2b2	Y132530
																															G2a2b2a1a1b1a1b1a2b2a	FT119236
																													G2a2b2a1a1b1a1b1a2c	FT204085
																														G2a2b2a1a1b1a1b1a2c1	FT245330
																										G2a2b2a1a1b1a1b2	BY27899.1
																											G2a2b2a1a1b1a1b2a	Z45474
																												G2a2b2a1a1b1a1b2a1	Z31348
																												G2a2b2a1a1b1a1b2a2	BY124520
																													G2a2b2a1a1b1a1b2a2a	Z45475
																													G2a2b2a1a1b1a1b2a2b	BY186574
																														G2a2b2a1a1b1a1b2a2b1	BY113713
																											G2a2b2a1a1b1a1b2b	BY113517
																												G2a2b2a1a1b1a1b2b1	BY55539
																										G2a2b2a1a1b1a1b3	BY180574
																									G2a2b2a1a1b1a1c	FGC809
																										G2a2b2a1a1b1a1c1	FGC807
																											G2a2b2a1a1b1a1c1a	FGC839
																												G2a2b2a1a1b1a1c1a2	BY108111
																											G2a2b2a1a1b1a1c1b	Z40867
																												G2a2b2a1a1b1a1c1b1	Z31846
																													G2a2b2a1a1b1a1c1b1a	Z31848
																														G2a2b2a1a1b1a1c1b1a2	Z43711
																													G2a2b2a1a1b1a1c1b1c	Y96077
																												G2a2b2a1a1b1a1c1b3	Z45501
																													G2a2b2a1a1b1a1c1b3a	Z45502
																														G2a2b2a1a1b1a1c1b3a1	BY117183
																												G2a2b2a1a1b1a1c1b5	BY11257
																													G2a2b2a1a1b1a1c1b5a	FT107767
																										G2a2b2a1a1b1a1c3	Z39863
																											G2a2b2a1a1b1a1c3a	BY128706
																												G2a2b2a1a1b1a1c3a1	Z45498
																													G2a2b2a1a1b1a1c3a1a	BY50750
																														G2a2b2a1a1b1a1c3a1a2	BY135965
																											G2a2b2a1a1b1a1c3b	BY60621
																												G2a2b2a1a1b1a1c3b1	BY200867
																								G2a2b2a1a1b1a2	Z6911
																									G2a2b2a1a1b1a2a	Y97776
																							G2a2b2a1a1b1b	S10458
																								G2a2b2a1a1b1b1	Z40133
																									G2a2b2a1a1b1b1a	Z24311
																										G2a2b2a1a1b1b1a1	Z24316
																										G2a2b2a1a1b1b1a2	S16124
																									G2a2b2a1a1b1b1c	BY152564
																						G2a2b2a1a1b2	Z27264
																							G2a2b2a1a1b2a	Z39088
																								G2a2b2a1a1b2a1	Z39525
																								G2a2b2a1a1b2a2	Z39531
																							G2a2b2a1a1b2b	Z39534
																					G2a2b2a1a1c	CTS342
																						G2a2b2a1a1c1	Z724
																							G2a2b2a1a1c1a	CTS4472
																								G2a2b2a1a1c1a1	CTS5990
																									G2a2b2a1a1c1a1a	CTS7045
																										G2a2b2a1a1c1a1a1	L640
																											G2a2b2a1a1c1a1a1a	Y88
																												G2a2b2a1a1c1a1a1a1	Z3307
																										G2a2b2a1a1c1a1a2	Z3428
																											G2a2b2a1a1c1a1a2a	Z26414
																												G2a2b2a1a1c1a1a2a1	FGC7477
																													G2a2b2a1a1c1a1a2a1a	Z27232
																														G2a2b2a1a1c1a1a2a1a1	Z6028
																															G2a2b2a1a1c1a1a2a1a1a	FGC249
																																G2a2b2a1a1c1a1a2a1a1a1	FGC263
																																G2a2b2a1a1c1a1a2a1a1a2	FGC31715
																																	G2a2b2a1a1c1a1a2a1a1a2a	FGC31721
																														G2a2b2a1a1c1a1a2a1a2	BY198799
																															G2a2b2a1a1c1a1a2a1a2a	BY200149
																												G2a2b2a1a1c1a1a2a2	FGC23437
																													G2a2b2a1a1c1a1a2a2a	FGC23440
																														G2a2b2a1a1c1a1a2a2a1	FGC23444
																														G2a2b2a1a1c1a1a2a2a2	Z36163
																															G2a2b2a1a1c1a1a2a2a2a	ALK539.2
																																G2a2b2a1a1c1a1a2a2a2a3	BY199033
																														G2a2b2a1a1c1a1a2a2a4	Z31300
																												G2a2b2a1a1c1a1a2b1	Z6434
																												G2a2b2a1a1c1a1a2b2	Z6033
																												G2a2b2a1a1c1a1a2b5	FT8419
																													G2a2b2a1a1c1a1a2b5a	FT99800
																											G2a2b2a1a1c1a1a2c	Z43085
																												G2a2b2a1a1c1a1a2c1	FGC58718
																													G2a2b2a1a1c1a1a2c1a	FGC58712
																										G2a2b2a1a1c1a1a4	BY79219
																										G2a2b2a1a1c1a1a5	BY49421
																										G2a2b2a1a1c1a1a7	FT304501
																									G2a2b2a1a1c1a1b	FGC46572
																										G2a2b2a1a1c1a1b1	FGC46568
																											G2a2b2a1a1c1a1b1a	FG46570
																												G2a2b2a1a1c1a1b1a2	Z44709
																						G2a2b2a1a1c2	FGC12126
																							G2a2b2a1a1c2a	L660
																								G2a2b2a1a1c2a1	Y40421
																								G2a2b2a1a1c2a2	FGC12134
																							G2a2b2a1a1c2b	Z16670
																								G2a2b2a1a1c2b1	Y165250
																									G2a2b2a1a1c2b1a	Z16713
																										G2a2b2a1a1c2b1a2	FGC64954
																											G2a2b2a1a1c2b1a2a	Z31214
																												G2a2b2a1a1c2b1a2a2	FT176549
																											G2a2b2a1a1c2b1a2c	FGC67660
																												G2a2b2a1a1c2b1a2c1	FGC67663
																													G2a2b2a1a1c2b1a2c1a	A20424.2
																									G2a2b2a1a1c2b1b	Z41276
																										G2a2b2a1a1c2b1b2	BY179507
																										G2a2b2a1a1c2b1b4	BY134284
																											G2a2b2a1a1c2b1b4a	FT165350
																										G2a2b2a1a1c2b1b5	FT184651
																								G2a2b2a1a1c2b2	FGC46781
																									G2a2b2a1a1c2b2a	Z31206
																										G2a2b2a1a1c2b2a1	FGC46784
																											G2a2b2a1a1c2b2a1c	Y109581
																										G2a2b2a1a1c2b2a2	Z45074
																								G2a2b2a1a1c2b4	FT232485
																						G2a2b2a1a1c3	PF4202
																							G2a2b2a1a1c3a	FGC37656
																								G2a2b2a1a1c3a1	FGC37724
																								G2a2b2a1a1c3a2	BY156830
																							G2a2b2a1a1c3b	BY144476
																						G2a2b2a1a1c4	FGC65058
																							G2a2b2a1a1c4b	BY21806
																					G2a2b2a1a1e	BY1433
																						G2a2b2a1a1e1	BY79491
																							G2a2b2a1a1e1a	BY63392
																						G2a2b2a1a1e2	Y81453
																							G2a2b2a1a1e2a	Z30708
																							G2a2b2a1a1e2b	FT273344
																					G2a2b2a1a1f	Z39263
																				G2a2b2a1a2	Z38302
																			G2a2b2a1b	Z30527
																				G2a2b2a1b1	Z31240
																				G2a2b2a1b2	PH1381
																		G2a2b2a2	M278
																		G2a2b2a3	Z6885
																			G2a2b2a3a	Z39310
																				G2a2b2a3a1	Z39308
																				G2a2b2a3a2	Z40745
																			G2a2b2a3b	Z6363
																				G2a2b2a3b1	Z6147
																			G2a2b2a4a	Z30519
																				G2a2b2a4a1	Z40458
																					G2a2b2a4a1a	Y21360
																						G2a2b2a4a1a1	Y21355
																							G2a2b2a4a1a1a	Z30522
																							G2a2b2a4a1a1b	BY92692
																								G2a2b2a4a1a1b2	CTS6435
																					G2a2b2a4a1b	BY56395
																						G2a2b2a4a1b2	BY62887
																							G2a2b2a4a1b2a	BY68738
																				G2a2b2a4a2	Z31387
																			G2a2b2a4b	BY96506
																	G2a2b2b	PF3359
																		G2a2b2b1	F1193
																			G2a2b2b1a	PF3369
																				G2a2b2b1a1	F872
																					G2a2b2b1a1a	PF3378
																						G2a2b2b1a1a1	PF3420
																										G2a2b2b1a1a1a1a1	Z45756
																						G2a2b2b1a1a2	Z7016
																							G2a2b2b1a1a2a	Z7022
																								G2a2b2b1a1a2a1	Z46195
																								G2a2b2b1a1a2a3	BY63384
																									G2a2b2b1a1a2a3b	BY77947
																								G2a2b2b1a1b1a1	M310.2
																								G2a2b2b1a1b1a2	FGC7263
																							G2a2b2b1a1b1b	BY4786. BY65351
																				G2a2b2b1a2	FGC82036
																			G2a2b2b1b	FGC52601
																				G2a2b2b1b1	BY193673
																				G2a2b2b1b2	FGC52622
																					G2a2b2b1b2c	FGC52629
																						G2a2b2b1b2c1	FGC52664
																		G2a2b2b2	PH488
																			G2a2b2b2a	CTS1455
													G2b	M3115
														G2b1	M377
															G2b1a	BY794
																G2b1a2	FGC32402
																	G2b1a2a	FGC35913
																		G2b1a2a1	FGC35915
																			G2b1a2a1a	Y67052
																		G2b1a2a2	Y16169
																			G2b1a2a2a	BY21363.2
																			G2b1a2a2b	BY165438
																		G2b1a2a3	BY147680
																		G2b1a2a4	BY37058
																	G2b1a2b	FGC32413
																		G2b1a2b1	FGC32409
																			G2b1a2b1a	FGC32397
																			G2b1a2b1b	FGC55702
																		G2b1a2b2	Y50032
																			G2b1a2b2b	BY181331
																	G2b1a2c	Z44684
														G2b2	FGC3022
															G2b2a	Z8022
															G2b2b	Z37343
											HIJK	F929
												H	L901
													H1	L902
														H1a	M69
															H1a1	M52
																H1a1a	M82
																	H1a1a1	Z14668
																		H1a1a1a	M36
																		H1a1a1b	L683
																	H1a1a4	M2914
																		H1a1a4a	Z14588
																		H1a1a4b	Z4361
																			H1a1a4b1	Z5873
																			H1a1a4b2	M2972
																				H1a1a4b2a	Z5876
																					H1a1a4b2a1	Z5877
																				H1a1a4b2b	Z5878
																					H1a1a4b2b1	Z5879
																				H1a1a4b2c	M3038
																					H1a1a4b2c1	Z5881
																						H1a1a4b2c1a	Z5882
																							H1a1a4b2c1a1	Z5883
																								H1a1a4b2c1a1a	Z5884
																			H1a1a4b3	Z4507
																				H1a1a4b3a	Z5886
																					H1a1a4b3a1	Z5887
																				H1a1a4b3b	Z5888
																					H1a1a4b3b1	Z5889
																						H1a1a4b3b1a	Z5890
																			H1a1a4b4	Z4489
																				H1a1a4b4a	Z4542
																H1a1b	Z4469
																	H1a1b1	Z4487
																		H1a1b1a	Z4417
																			H1a1b1a1	Z14686
																		H1a1b1b	Z40935
															H1a2	Z5867
																H1a2a	Apt
																H1a2b	Z14258
																	H1a2b1	Z14308
																		H1a2b1a	Z5868
																		H1a2b1b	Z34492
																		H1a2b1c	B368
													H2	P96
														H2a	Z19067
																H2a1a	Z19049
															H2b1	Y21600
													H3	Z5857
														H3a	Z5866
															H3a1	Z5864
															H3a2	Z5863
																H3a2a	Z5865
																	H3a2a1	Z12646
																		H3a2a1a	Y11965
																	H3a2a2	Z5862
																H3a2b	Z5858
														H3b	Z13871
															H3b1	Z5859
															H3b2	PH24
												IJK	L15
													IJ	M429
														I	M170
															I1	M253
																I1a	DF29
																	I1a1	CTS6364
																						I1a1a1a1a	M227
																									I1a1a1a1a2a1	M72
																		I1a1b	CTS10028
																			I1a1b1	L22
																					I1a1b1a1	P109
																					I1a1b1a3	L205.1
																						I1a1b1a4a	Z74
																								I1a1b1a4a1a	L287
																									I1a1b1a4a1a1	L258
																										I1a1b1a4a1a1b	CTS2242
																												I1a1b1a4a1a1b1a	Z134
																										I1a1b1a4a1a1c	Z2045
																							I1a1b1a4a2	L813
																							I1a1b1b1a1	L300
																			I1a1b2	FGC10477
																	I1a2	S244
																		I1a2a	S246
																				I1a2a1a	Z62
																						I1a2a1a1a	S440
																							I1a2a1a1a1	S1953
																								I1a2a1a1a1a	S1954
																									I1a2a1a1a1a1	L338
																								I1a2a1a1a1b	CTS10937
																							I1a2a1a1a2	F2642.1
																								I1a2a1a1a2a	CTS6772
																							I1a2a1a1a3	A196
																								I1a2a1a1a3a	Y6900
																									I1a2a1a1a3a1	Y7140
																										I1a2a1a1a3a1a	Y7477
																										I1a2a1a1a3a2a	Y6885
																						I1a2a1a1d	CTS7362
																								I1a2a1a1d1a	S247
																									I1a2a1a1d1a1	BY266
																										I1a2a1a1d1a1a	L1302
																											I1a2a1a1d1a1a1	BY147
																											I1a2a1a1d1a1a2	BY126
																								I1a2a1a1d1b	L573
																							I1a2a1a1d2	L1248
																										I1a2a1a1d2a1c	L803
																			I1a2a2	Z2041
																				I1a2a2a	Z2040
																					I1a2a2a3	Y2170
																						I1a2a2a3a	Z2042
																		I1a2b	S296.1
																			I1a2b1	Z2541
																	I1a3	S243
																		I1a3a	S2078
																				I1a3a1a	Y2245.2
																					I1a3a1a1	L1237
																						I1a3a1a1a	FGC9550
																					I1a3a1a2	S10360
																						I1a3a1a2a	S15301
																							I1a3a1a2a1	Y6228
																					I1a3a1a3	FGC14479
																				I1a3a1b	Y6375
																		I1a3b	BY351
																			I1a3b1	BY332
																				I1a3b1a	Y7060
																					I1a3b1a1	BY11
																				I1a3b1b	Y13946
																			I1a3b2	CTS10345
																					I1a3b2a1	Y10994
																			I1a3b3	CTS4279
																		I1a3c	BY62
																I1b	S249
																	I1b1	CTS6397
																I1c	Y18119
															I2	M438
																I2a	CTS1799
																	I2a1	L460
																		I2a1a	P37.2
																			I2a1a1	CTS595
																				I2a1a1a	M26
																							I2a1a1a1a1	L160
																								I2a1a1a1a1a	PF4088
																										I2a1a1a1a1a1a	CTS11338
																											I2a1a1a1a1a1a1	PF4189
																											I2a1a1a1a1a1a2	Z118
																								I2a1a1a1a1b	F1295
																				I2a1a1b	S21825
																					I2a1a1b1	L1286
																						I2a1a1b1a	L1287
																							I2a1a1b1a1	L233
																									I2a1a1b1a1a1	A417
																						I2a1a1b1b	L880
																					I2a1a1b2	L1294
																			I2a1a2	M423
																				I2a1a2a	L161.1
																					I2a1a2a1	AMM047
																						I2a1a2a1a	AMM079.1
																							I2a1a2a1a1	AMM068
																								I2a1a2a1a1a	AMM083
																									I2a1a2a1a1a1	A16983
																									I2a1a2a1a1a2	A2330
																										I2a1a2a1a1a2a	A2293
																											I2a1a2a1a1a2a1	CTS8849
																									I2a1a2a1a1a3	CTS4122
																									I2a1a2a1a1a4	A1513
																										I2a1a2a1a1a4a	A17998
																											I2a1a2a1a1a4a1	A16449
																								I2a1a2a1a1b	FGC7206
																									I2a1a2a1a1b1	S7708
																									I2a1a2a1a1b2	A10512
																										I2a1a2a1a1b2a	A11113
																									I2a1a2a1a1b3	FGC7184
																										I2a1a2a1a1b3a	FGC14451
																											I2a1a2a1a1b3a1	A7729
																												I2a1a2a1a1b3a1a	A12374
																											I2a1a2a1a1b3a2	FGC14457
																												I2a1a2a1a1b3a2a	FGC14463
																										I2a1a2a1a1b3b	FGC7156
																											I2a1a2a1a1b3b1	FGC14714
																												I2a1a2a1a1b3b1a	FGC14715
																											I2a1a2a1a1b3b2	FGC7218
																												I2a1a2a1a1b3b2a	A10741
																												I2a1a2a1a1b3b2b	A14890
																							I2a1a2a1a2	A1150
																								I2a1a2a1a2a	A8742
																									I2a1a2a1a2a1	A18001
																								I2a1a2a1a2b	A12014
																									I2a1a2a1a2b1	A11114
																									I2a1a2a1a2b2	A10036
																										I2a1a2a1a2b2a	A18004
																							I2a1a2a1a3	A10028
																								I2a1a2a1a3a	A14889
																						I2a1a2a1b	Y12072
																							I2a1a2a1b1	PF4135
																								I2a1a2a1b1a	Y11772
																									I2a1a2a1b1a1	A11115
																									I2a1a2a1b1a2	Y12075
																										I2a1a2a1b1a2a	A8270.2
																										I2a1a2a1b1a2b	Y12059
																											I2a1a2a1b1a2b1	Y19285
																											I2a1a2a1b1a2b2	A11984
																								I2a1a2a1b1b	A11374
																									I2a1a2a1b1b1	A17708
																								I2a1a2a1b1c	A13665.2
																									I2a1a2a1b1c1	A13664
																										I2a1a2a1b1c1a	A13904
																					I2a1a2a2	Y13338
																						I2a1a2a2a	Y13335
																							I2a1a2a2a1	Y13568
																								I2a1a2a2a1a	Y13589
																				I2a1a2b	L621
																					I2a1a2b1	CTS10936
																						I2a1a2b1a	CTS4002
																							I2a1a2b1a1	CTS5966
																								I2a1a2b1a1a	S9952
																									I2a1a2b1a1a1	S17250
																										I2a1a2b1a1a1a	BY128
																											I2a1a2b1a1a1a1	Z16971
																												I2a1a2b1a1a1a1a	A815
																													I2a1a2b1a1a1a1a1	A5875
																														I2a1a2b1a1a1a1a1a	A5874
																															I2a1a2b1a1a1a1a1a1	A5876
																														I2a1a2b1a1a1a1a1b	A14798
																												I2a1a2b1a1a1a1b	A2423
																												I2a1a2b1a1a1a1c	A16681
																											I2a1a2b1a1a1a2	A14973
																										I2a1a2b1a1a1b	Y4882
																											I2a1a2b1a1a1b2	B474
																											I2a1a2b1a1a1b3	Y15928
																												I2a1a2b1a1a1b3a	A7318
																											I2a1a2b1a1a1b4	A12505
																										I2a1a2b1a1a1c	PH908
																											I2a1a2b1a1a1c1	A356
																												I2a1a2b1a1a1c1a	A493
																													I2a1a2b1a1a1c1a1	Y6651
																													I2a1a2b1a1a1c1a2	A8740
																											I2a1a2b1a1a1c2	A5913
																												I2a1a2b1a1a1c2a	A16863
																											I2a1a2b1a1a1c3	PH3310
																											I2a1a2b1a1a1c4	A13912
																									I2a1a2b1a1a2	Y4460
																										I2a1a2b1a1a2a	Y3118
																											I2a1a2b1a1a2a1	Y5598
																												I2a1a2b1a1a2a1a	CTS5779
																													I2a1a2b1a1a2a1a1	Y10622
																										I2a1a2b1a1a2b	S8201
																											I2a1a2b1a1a2b1	Y13498
																										I2a1a2b1a1a2c	A6105
																											I2a1a2b1a1a2c1	Y16810
																									I2a1a2b1a1a3	Z17855
																										I2a1a2b1a1a3a	A1221
																										I2a1a2b1a1a3b	A16413
																									I2a1a2b1a1a4	Y18331
																										I2a1a2b1a1a4a	A7134
																											I2a1a2b1a1a4a1	A14877
																										I2a1a2b1a1a4b	A10959
																							I2a1a2b1a2	FGC20479
																								I2a1a2b1a2a	A16633
																								I2a1a2b1a2b	FGC20473
																					I2a1a2b2	A17060
																		I2a1b	M436
																			I2a1b1	M223
																				I2a1b1a	CTS616
																					I2a1b1a1	FGC15073
																						I2a1b1a1a	M284
																							I2a1b1a1a1	L1195
																								I2a1b1a1a1a	L126
																									I2a1b1a1a1a1	FGC20048
																											I2a1b1a1a1a1a1	FGC20065
																													I2a1b1a1a1a1a1a1	Y10655
																												I2a1b1a1a1a1a1b	FGC20064
																													I2a1b1a1a1a1a1b1	FGC20078
																														I2a1b1a1a1a1a1b1b	Y5673
																															I2a1b1a1a1a1a1b1b1	FGC21571
																														I2a1b1a1a1a1a1b1c	Y8599
																											I2a1b1a1a1a1a2	Y7190
																								I2a1b1a1a1b	L1193
																										I2a1b1a1a1b1a	CTS4922
																						I2a1b1a1b	Z2057
																							I2a1b1a1b1	L1229
																								I2a1b1a1b1a	Y3681
																										I2a1b1a1b1a1a	Z2054
																											I2a1b1a1b1a1a1	BY524
																												I2a1b1a1b1a1a1a	L812
																													I2a1b1a1b1a1a1a1	Y5308
																													I2a1b1a1b1a1a1a2	Y13008
																												I2a1b1a1b1a1a1b	Y10648
																											I2a1b1a1b1a1a2	Y7244
																												I2a1b1a1b1a1a2a	P53.3
																												I2a1b1a1b1a1a2b	Y7243
																											I2a1b1a1b1a1a3	FGC15106
																													I2a1b1a1b1a1a3a1	Y4760
																											I2a1b1a1b1a1a4	BY138
																								I2a1b1a1b1b	S18331
																									I2a1b1a1b1b1	L1230
																							I2a1b1a1b2	Y7240
																					I2a1b1a2	CTS10057
																						I2a1b1a2a	L701
																							I2a1b1a2a1	P78
																								I2a1b1a2a1a	S25733
																									I2a1b1a2a1a1	A427
																										I2a1b1a2a1a1a	S23612
																													I2a1b1a2a1a1a1a1	Y5360
																														I2a1b1a2a1a1a1a1a	Y5369
																											I2a1b1a2a1a1a2	S10702
																												I2a1b1a2a1a1a2a	Y11261
																								I2a1b1a2a1b	Y7219
																									I2a1b1a2a1b1	Y7214
																									I2a1b1a2a1b2	Y8945
																								I2a1b1a2a2a	L699
																									I2a1b1a2a2a1	L704
																											I2a1b1a2a2a1a1	Y7641
																									I2a1b1a2a2a2	S12195
																											I2a1b1a2a2a2a1	Y6973
																												I2a1b1a2a2a2a1a	BY159
																						I2a1b1a2b	Z161
																							I2a1b1a2b1	L801
																									I2a1b1a2b1a1	CTS1977
																										I2a1b1a2b1a1a	FGC33292
																												I2a1b1a2b1a1a1a	S8522
																													I2a1b1a2b1a1a1a1	P95
																											I2a1b1a2b1a1a2	CTS1858
																												I2a1b1a2b1a1a2a	CTS10148
																												I2a1b1a2b1a1b1a	S16196
																													I2a1b1a2b1a1b1a1	Y7152
																										I2a1b1a2b1a1c	BY526
																									I2a1b1a2b1a2	CTS6433
																										I2a1b1a2b1a2a	FGC3622
																											I2a1b1a2b1a2a1	FGC3618
																												I2a1b1a2b1a2a1a	Z78
																													I2a1b1a2b1a2a1a1	CTS8584
																														I2a1b1a2b1a2a1a1a	FGC3621
																															I2a1b1a2b1a2a1a1a1	FGC3620
																																I2a1b1a2b1a2a1a1a1a	L1198
																																	I2a1b1a2b1a2a1a1a1a1	S20905
																																		I2a1b1a2b1a2a1a1a1a1a	Z190
																																			I2a1b1a2b1a2a1a1a1a1a1	S434
																																				I2a1b1a2b1a2a1a1a1a1a1a	Y5729
																																				I2a1b1a2b1a2a1a1a1a1a1b	PH3769
																																				I2a1b1a2b1a2a1a1a1a1a1c	Y8712
																																				I2a1b1a2b1a2a1a1a1a1a2a	Y7280
																																	I2a1b1a2b1a2a1a1a1a2	P195.2
																																	I2a1b1a2b1a2a1a1a1a3	Y6060
																																		I2a1b1a2b1a2a1a1a1a3a	Y5748
																																			I2a1b1a2b1a2a1a1a1a3a1	Y7272
																																		I2a1b1a2b1a2a1a1a1a3b	Y7273
																												I2a1b1a2b1a2a1b	Y9161
																											I2a1b1a2b1a2a2	Y4925
																												I2a1b1a2b1a2a2a	CTS661
																													I2a1b1a2b1a2a2a1	Y5717
																												I2a1b1a2b1a2a2b	FGC17399
																														I2a1b1a2b1a2a2b1a	S8104
																											I2a1b1a2b1a2a3	FGC19998
																												I2a1b1a2b1a2a3a	Y5695
																													I2a1b1a2b1a2a3a1	Y7263
																														I2a1b1a2b1a2a3a1a	Y7265
																													I2a1b1a2b1a2a3a2	FGC20004
																														I2a1b1a2b1a2a3b1a	Y10659
																															I2a1b1a2b1a2a3b1a1	Y11231
																										I2a1b1a2b1a2b	S25383
																												I2a1b1a2b1a2b1a	CTS5332.2
																												I2a1b1a2b1a2b2a	Y4769
																										I2a1b1a2b1a2c	Y7426
																										I2a1b1a2b1a3a	L1290
																									I2a1b1a2b1a4	Y7202
																							I2a1b1a2b2	L623
																							I2a1b1a2b3	CTS11871
																								I2a1b1a2b3a	Y10671
																				I2a1b1b	S9403
																						I2a1b1b1a	L1228
																					I2a1b1b2	Y6099
																			I2a1b2	FGC29562
																				I2a1b2a	L38
																					I2a1b2a1	BY14072
																						I2a1b2a1a	L533
																							I2a1b2a1a1	Y52464
																						I2a1b2a1b	CTS10544
																					I2a1b2a2	S2606
																						I2a1b2a2a	S24121
																							I2a1b2a2a1	S19763
																								I2a1b2a2a1a	Y16415
																									I2a1b2a2a1a1	Y17121
																										I2a1b2a2a1a1a	Y17117
																											I2a1b2a2a1a1a1	BY14018
																												I2a1b2a2a1a1a1a	Y17266
																										I2a1b2a2a1a1b	BY14650
																											I2a1b2a2a1a1b1	BY14649
																											I2a1b2a2a1a1b2	Y86214
																												I2a1b2a2a1a1b2a	Y89787
																									I2a1b2a2a1a2	Y16417
																								I2a1b2a2a1b	FGC36959
																									I2a1b2a2a1b1	BY14023
																										I2a1b2a2a1b1a	FGC36962
																											I2a1b2a2a1b1a1	BY25360
																												I2a1b2a2a1b1a1a	BY25379
																											I2a1b2a2a1b1a2	FGC36963
																							I2a1b2a2a2	Y86523
																								I2a1b2a2a2a	BY42260
																								I2a1b2a2a2b	BY44566
																							I2a1b2a2a3	BY37418
																								I2a1b2a2a3a	BY44939
																						I2a1b2a2b	PH1237
																							I2a1b2a2b1	PH2591
																								I2a1b2a2b1a	Y32631
																									I2a1b2a2b1a1	Y33743
																							I2a1b2a2b2	BY14026
																								I2a1b2a2b2a	BY14037
																									I2a1b2a2b2a1	BY25362
																								I2a1b2a2b2b	BY25359
																									I2a1b2a2b2b1	A577
																						I2a1b2a2c	S2488
																							I2a1b2a2c1	BY14048
																							I2a1b2a2c2	FGC29631
																								I2a1b2a2c2a	S8239
																									I2a1b2a2c2a1	S19023
																						I2a1b2a2d	BY1183
																							I2a1b2a2d1	Y18919
																								I2a1b2a2d1a	S4556
																									I2a1b2a2d1a1	Y20294
																									I2a1b2a2d1a2	BY14000
																										I2a1b2a2d1a2a	BY14011
																									I2a1b2a2d1a3	Y105542
																								I2a1b2a2d1b	Y18921
																									I2a1b2a2d1b1	BY4147
																										I2a1b2a2d1b1a	BY4146
																										I2a1b2a2d1b1b	Y20202
																	I2a2	L596
																		I2a2a	PF6915
																				I2a2a1a	L1251.1
																					I2a2a1a1	S14736
																						I2a2a1a1a	FGC18545
																							I2a2a1a1a1	FGC18649
																							I2a2a1a1a2	F313
																						I2a2a1a1b	S18142
																				I2a2a1b	S6596
																					I2a2a1b1	S6595
																						I2a2a1b1a	BY2801
																							I2a2a1b1a1	A1144
																						I2a2a1b1b	CTS4092
																							I2a2a1b1b1	FGC18098
																								I2a2a1b1b1a	BY23
																									I2a2a1b1b1a1	BY19
																										I2a2a1b1b1a1a	BY1265
																											I2a2a1b1b1a1a1	F2044
																												I2a2a1b1b1a1a1a	BY43.3
																										I2a2a1b1b1a1b	BY48
																								I2a2a1b1b1b	BY2797
																									I2a2a1b1b1b1	BY25
																										I2a2a1b1b1b1a	BY24
																								I2a2a1b1b1c	L199.2
																					I2a2a1b2	BY45868
																			I2a2a2	BY431
																		I2a2b	PH2569
																			I2a2b1	Z26403
																					I2a2b1b2	BY16408
																				I2a2b1d	BY2808
																					I2a2b1d1	BY2807
																I2b	L415
																	I2b1	BY32138
														J	M304
															J1	L255
																J1a	CTS5368
																				J1a2a1a	FGC14316
																						J1a2a1a1a	P56
																					J1a2a1a2	P58
																						J1a2a1a2a	L92.1
																						J1a2a1a2b	L147.1
																							J1a2a1a2c1	L817
																								J1a2a1a2c1a	L818
																									J1a2a1a2c1a1	L816
																																			J1a2a1a2d2b2b2c4d2a2a5	L222.2
														"	"
																		J1a2b	CTS15
																			J1a2b1	Z1842
																J1b	F1614
															J2	M172
																J2a	M410
																	J2a1	PF4610
																		J2a1a	F4326
																			J2a1a1	PF5088
																				J2a1a1a	PF5125
																					J2a1a1a2	Z2229
																						J2a1a1a2a	Z6065
																							J2a1a1a2a2	Y8522
																								J2a1a1a2a2a	Z39478
																									J2a1a1a2a2a2	P81
																								J2a1a1a2a2b	M47
																						J2a1a1a2b	CTS4800
																							J2a1a1a2b1	S23154
																								J2a1a1a2b1b	M319
																							J2a1a1a2b2	M67
																								J2a1a1a2b2a	PF5126
																										J2a1a1a2b2a1a	M92
																												J2a1a1a2b2a1a1a	S16400
																													J2a1a1a2b2a1a1a3	Y20051
																														J2a1a1a2b2a1a1a3b	L560
																									J2a1a1a2b2a2	FGC7860
																										J2a1a1a2b2a2b	FGC7861
																											J2a1a1a2b2a2b1	FGC45722
																													J2a1a1a2b2a2b1b2	M166
																											J2a1a1a2b2a2b3	FGC3291
																												J2a1a1a2b2a2b3a	L210
																				J2a1a1b	AM01359
																					J2a1a1b1	PF5197
																						J2a1a1b1a	PF5172
																							J2a1a1b1a1	Z7314
																								J2a1a1b1a1a	S21554
																									J2a1a1b1a1a1	FGC16141
																										J2a1a1b1a1a1a	S15439
																											J2a1a1b1a1a1a1	S24735
																										J2a1a1b1a1a1b	FGC16143
																					J2a1a1b2	F4168
																						J2a1a1b2a	AM00103
																							J2a1a1b2a1	L25
																								J2a1a1b2a1a	F3133.1
																									J2a1a1b2a1a1	A5360
																										J2a1a1b2a1a1a	M158.2
																									J2a1a1b2a1a2	FGC9878
																										J2a1a1b2a1a2a	SK1382
																											J2a1a1b2a1a2a1	FGC30637.2
																												J2a1a1b2a1a2a1a	L271
																								J2a1a1b2a1b	FGC2837
																									J2a1a1b2a1b1	L70
																													J2a1a1b2a1b1b2b1	M318
																									J2a1a1b2a1c1	L243
																										J2a1a1b2a1c2a	L254
																									J2a1a2b2a2a2	M68
																		J2a2a	L581
																							J2a2a1a1a2	P279
																									J2a2a1a1a2b1	M340
																J2b	M12
																	J2b1	M205
																		J2b2a	M241
																			J2b2a1	L283
																									J2b2a1a1a1a1	AMM450
																										J2b2a1a1a1a1a	Z1297
																												J2b2a1a1a1a1a1a	Z631
																					J2b2a2b2	Z2437
																						J2b2a2b2a	Y965
																							J2b2a2b2a1	Z8347
																						J2b2a2b2b	Z2443
																							J2b2a2b2b1	Z2448
													K	M9
														LT	L298
															L	M20
																L1	M22
																	L1a	M2481
																		L1a1	M27
																			L1a1a	Z20458
																				L1a1a2	Z20475
																			L1a1b	Z5926
																				L1a1b1	Z20519
																				L1a1b2	Z20520
																					L1a1b2a	Z5929
																				L1a1b3	Z5930
																					L1a1b3a	Z5931
																						L1a1b3a1	Z5933
																							L1a1b3a1a	Z5934
																								L1a1b3a1a1	Z20500
																								L1a1b3a1a2	Z8033
																									L1a1b3a1a2a	Z5937
																										L1a1b3a1a2a1	Z5938
																						L1a1b3a2	Z20568
																		L1a2	M357
																			L1a2a	M2398
																				L1a2a1	Z5920
																					L1a2a1a	Z20284
																						L1a2a1a1	Z8030
																							L1a2a1a1a	Z5923
																	L1b	FGC36877
																		L1b1	M349
																			L1b1a	Page116
																		L1b2	SK1412
																			L1b2a	M274
																			L1b2b	PH8
																				L1b2b1	PH2079
																					L1b2b1a	PH1714
																				L1b2b2	PH1099
																					L1b2b2a	PH1728
																					L1b2b2b	Y18889
																			L1b2c	SK1414
																		L1b3	CTS7154.2
																			L1b3a	CTS3080
																L2	L595
															T	M184
																T1	L206
																	T1a	M70
																		T1a1	L162
																			T1a1a	L208
																				T1a1a1	CTS11451
																					T1a1a1a	PF7455
																						T1a1a1a1	Y12643
																							T1a1a1a1a	Y12641
																								T1a1a1a1a1	Y15122
																							T1a1a1a1b	PF7433
																								T1a1a1a1b1	PF7445
																									T1a1a1a1b1a	PF7444
																					T1a1a1b	FGC3995
																						T1a1a1b1	Y16244
																							T1a1a1b1a	FGC40894
																								T1a1a1b1a1	Y14121
																						T1a1a1b2	CTS2214
																							T1a1a1b2a	Y15127
																								T1a1a1b2a1	Y15711
																									T1a1a1b2a1a	Y15709
																									T1a1a1b2a1b	FGC41002
																								T1a1a1b2a2	Y3781
																									T1a1a1b2a2a	Y3824
																										T1a1a1b2a2a1	Y3782
																											T1a1a1b2a2a1a	Y3836
																										T1a1a1b2a2a2	Y9102
																											T1a1a1b2a2a2a	Y16247
																											T1a1a1b2a2a2b	Y9104
																												T1a1a1b2a2a2b1	Y16241
																													T1a1a1b2a2a2b1a	Y16240
																												T1a1a1b2a2a2b2	FGC40425
																													T1a1a1b2a2a2b2a	FGC40529
																							T1a1a1b2b	Z709
																								T1a1a1b2b1	Y9109
																									T1a1a1b2b1a	Y18004
																										T1a1a1b2b1a1	Y18003
																								T1a1a1b2b2	FGC3988
																									T1a1a1b2b2a	Y6409
																									T1a1a1b2b2b	CTS660
																										T1a1a1b2b2b1	CTS8603
																											T1a1a1b2b2b1a	P77
																												T1a1a1b2b2b1a1	CTS6174
																													T1a1a1b2b2b1a1a	CTS6507
																														T1a1a1b2b2b1a1a1	CTS9882
																															T1a1a1b2b2b1a1a1a	Y15648
																																T1a1a1b2b2b1a1a1a1	FGC40667
																															T1a1a1b2b2b1a1a1b	FGC4018
																																T1a1a1b2b2b1a1a1b1	Y16136
																																T1a1a1b2b2b1a1a1b2	FGC4040
																																	T1a1a1b2b2b1a1a1b2a	FGC4065
																																		T1a1a1b2b2b1a1a1b2a1	FGC4016
																															T1a1a1b2b2b1a1a1c	CTS6901
																														T1a1a1b2b2b1a1a2	CTS6280
																														T1a1a1b2b2b1a1a3	Y16021
																														T1a1a1b2b2b1a1a4	Y11696
																				T1a1a2	Y16897
																					T1a1a2a	Z19963
																		T1a2	L131
																			T1a2a	PH141
																				T1a2a1	P322
																			T1a2b	L446
																				T1a2b1	CTS933
																					T1a2b1a	CTS3767
																						T1a2b1a1	CTS8862
																							T1a2b1a1a	Y17493
																				T1a2b2	FGC23010
																					T1a2b2a	FGC23024
																		T1a3	FGC1350
																			T1a3a	Y11675
																			T1a3b	FGC1340
																				T1a3b1	L1255
																					T1a3b1a	FGC1408
																				T1a3b2	Y13280
																					T1a3b2a	Y13698
																						T1a3b2a1	SK1474
																							T1a3b2a1a	FGC41095
																							T1a3b2a1b	Y13293
																								T1a3b2a1b1	Y13279
																T2	PH110
														NO	F549
															NO1	M214
																N	"M231
																CTS6605/M2214, CTS7230/M2221,CTS7264/M2222, CTS7433/M2224, CTS7551/F3708/L667/M2226,  CTS8687/M2241, CTS8893, CTS9197/F3714/M2246, CTS9920/M2251, CTS10527/M2258, CTS10559/M2259,  CTS11117/M2284/Z5028, CTS11321^/M2287^, CTS11491/M2290, F1052/M2144, F1053/M2145/S10579/V1392, F1359/M2158/S11735/V2374, F1715/M2168/S15571, F1815/M2177/V2528, F1951/M2189/V2951, F2080/M2199/S17957, F2088/M2200,  F2201/M2205/S18745, F2344^/M2211^, F2509/M2225/S19819, F2571/M2230/S20099, F2596^/M2233^/S20231^, F2636/L734/M2237, F2692^/M2239^/S20571^, F2783/M2245, F2919/M2250/V3683, F2968, F2981/M2252/S21927/V3903, F2999/M2253/S22032/V3957, F3108/M2260, F3123/M2262/S22586,  F3227/M2100/S23418, F3235^/M2270^, F3299/M2275, F3308/M2277/S24122, F3422/M2294, F3426/M2296/S25467, F3691/M2148/V1495, F3693/M2160/V2508, F3694/M2161, M2151/S11174/V1731/Z4852, M2159/S11864/V2504/Z4861, M2263/Z4931, M2279/S24429/Z5023, M2280/Z5025, M2281/Z5026 "	CTS6605/M2214, CTS7230/M2221,CTS7264/M2222, CTS7433/M2224, CTS7551/F3708/L667/M2226,  CTS8687/M2241, CTS8893, CTS9197/F3714/M2246, CTS9920/M2251, CTS10527/M2258, CTS10559/M2259,  CTS11117/M2284/Z5028, CTS11321^/M2287^, CTS11491/M2290, F1052/M2144, F1053/M2145/S10579/V1392, F1359/M2158/S11735/V2374, F1715/M2168/S15571, F1815/M2177/V2528, F1951/M2189/V2951, F2080/M2199/S17957, F2088/M2200,  F2201/M2205/S18745, F2344^/M2211^, F2509/M2225/S19819, F2571/M2230/S20099, F2596^/M2233^/S20231^, F2636/L734/M2237, F2692^/M2239^/S20571^, F2783/M2245, F2919/M2250/V3683, F2968, F2981/M2252/S21927/V3903, F2999/M2253/S22032/V3957, F3108/M2260, F3123/M2262/S22586,  F3227/M2100/S23418, F3235^/M2270^, F3299/M2275, F3308/M2277/S24122, F3422/M2294, F3426/M2296/S25467, F3691/M2148/V1495, F3693/M2160/V2508, F3694/M2161, M2151/S11174/V1731/Z4852, M2159/S11864/V2504/Z4861, M2263/Z4931, M2279/S24429/Z5023, M2280/Z5025, M2281/Z5026 "
																	N1	CTS3750
																		N1a	F1206
																			N1a1	M46
																				N1a1a	M178
																						N1a1a1a	L708
																							N1a1a1a1	P298
																								N1a1a1a1a	L392
																									N1a1a1a1a1	CTS10760
																										N1a1a1a1a1a	CTS2929
																												N1a1a1a1a1a1a	L550
																													N1a1a1a1a1a1a1	B215
																														N1a1a1a1a1a1a1a	M2784.1
																															N1a1a1a1a1a1a1a1	L551
																															N1a1a1a1a1a1a1a2	BY158
																																N1a1a1a1a1a1a1a2a	L591
																																N1a1a1a1a1a1a1a2b	L1027
																															N1a1a1a1a1a1a1a3	FGC13372
																															N1a1a1a1a1a1a1a4	CTS8173
																														N1a1a1a1a1a1a1b	Y4706
																													N1a1a1a1a1a1a3	Y9455
																											N1a1a1a1a1a2	CTS9976
																												N1a1a1a1a1a2a	Y5002
																														N1a1a1a1a1a2a1a	CTS2601
																															N1a1a1a1a1a2a1a1	PH3744
																																	N1a1a1a1a1a2a1a1a2	PH547
																															N1a1a1a1a1a2a1a2	FGC14180
																									N1a1a1a1a2	Z1936
																											N1a1a1a1a2a1	CTS2733
																												N1a1a1a1a2a1a	CTS10035
																														N1a1a1a1a2a1a1a	CTS1737
																															N1a1a1a1a2a1a1a1	Z19813
																																N1a1a1a1a2a1a1a1a	CTS1950
																																	N1a1a1a1a2a1a1a1a1	Z1941
																																		N1a1a1a1a2a1a1a1a1a	Z1940
																																			N1a1a1a1a2a1a1a1a1a1	Z5033
																																				N1a1a1a1a2a1a1a1a1a1a	Z5038
																																					N1a1a1a1a2a1a1a1a1a1a1	CTS4303
																																						N1a1a1a1a2a1a1a1a1a1a1a	CTS8445
																																							N1a1a1a1a2a1a1a1a1a1a1a1	CTS12908
																																N1a1a1a1a2a1a1a1b	CTS1152
																																		N1a1a1a1a2a1a1a1b1a	Z5892
																																			N1a1a1a1a2a1a1a1b1b1	Z5893
																																				N1a1a1a1a2a1a1a1b1b1a	YP1141
																									N1a1a1a1a3	B197
																										N1a1a1a1a3a	F4205
																											N1a1a1a1a3a2	B219
																												N1a1a1a1a3a2c	Z35327
																													N1a1a1a1a3a2c2	B199
																										N1a1a1a1a3b	B202
																											N1a1a1a1a3b1	B203
																											N1a1a1a1a3b2	B204
																												N1a1a1a1a3b2b	B206
																													N1a1a1a1a3b2b2	B208
																									N1a1a1a1a4	M2019
																											N1a1a1a1a4a1	M1982
																												N1a1a1a1a4a1a	M1979
																													N1a1a1a1a4a1a1	M1984
																							N1a1a1a2	B211
																								N1a1a1a2b	B181
																									N1a1a1a2b2	B229
																			N1a2	F1008
																				N1a2a	M128
																				N1a2b	B523
																					N1a2b1	P63
																						N1a2b1b	B169
																							N1a2b1b1	B170
																								N1a2b1b1b	B172
																							N1a2b1b2	B175
																								N1a2b1b2a	B227
																									N1a2b1b2a2	Z35131
																								N1a2b1b2b	PF3415
																						N1a2b2a	FGC10821.2
																							N1a2b2a1	VL97
																						N1a2b2b	Y23786
																					N1a2b3	B525
																		N1b	F2930
																			N1b1	CTS582
																				N1b1a	Y6374
																					N1b1a1	CTS7324
																			N1b2	M1819
																				N1b2a	M1811
																O	M175
																	O1	F265
																		O1a	M119
																			O1a1	B384
																				O1a1a	M307.1
																					O1a1a1	F446
																						O1a1a1a	F140
																							O1a1a1a1	F78
																								O1a1a1a1a	F81
																									O1a1a1a1a1	CTS2458
																										O1a1a1a1a1a	F533
																											O1a1a1a1a1a1	F492
																												O1a1a1a1a1a1a	F656
																													O1a1a1a1a1a1a1	A12440
																														O1a1a1a1a1a1a1a	A12439
																															O1a1a1a1a1a1a1a1	MF20192
																															O1a1a1a1a1a1a1a2	MF6492
																														O1a1a1a1a1a1a1b	MF6062
																														O1a1a1a1a1a1a1c	MF6067
																															O1a1a1a1a1a1a1c1	MF6065
																																O1a1a1a1a1a1a1c1a	MF6277
																														O1a1a1a1a1a1a1d	MF6409
																														O1a1a1a1a1a1a1e	A23676. MF20133
																													O1a1a1a1a1a1a2	A14788
																													O1a1a1a1a1a1a3	F65
																														O1a1a1a1a1a1a3a	Y137223
																													O1a1a1a1a1a1a4	MF1069
																														O1a1a1a1a1a1a4a	MF1068
																													O1a1a1a1a1a1a5	Z23482
																													O1a1a1a1a1a1a6	MF2651
																													O1a1a1a1a1a1a7	Y137090
																														O1a1a1a1a1a1a7a	Y137081
																													O1a1a1a1a1a1a8	MF6069
																													O1a1a1a1a1a1a9	MF6073
																													O1a1a1a1a1a1a10	F18460
																														O1a1a1a1a1a1a10a	F19365
																													O1a1a1a1a1a1a11	ACT1077
																													O1a1a1a1a1a1a12	ACT3892
																													O1a1a1a1a1a1a13	ACT6705
																														O1a1a1a1a1a1a13a	ACT6493
																												O1a1a1a1a1a1b	FGC66168
																													O1a1a1a1a1a1b1	CTS11553
																														O1a1a1a1a1a1b1a	MF14611
																															O1a1a1a1a1a1b1a1	MF14730
																															O1a1a1a1a1a1b1a2	MF14734
																														O1a1a1a1a1a1b1b	MF15219
																															O1a1a1a1a1a1b1b1	MF16753
																																O1a1a1a1a1a1b1b1a	MF15176
																																O1a1a1a1a1a1b1b1b	MF14234
																														O1a1a1a1a1a1b1c	FGC66158
																															O1a1a1a1a1a1b1c1	FGC66159
																																O1a1a1a1a1a1b1c1a	MF15014
																															O1a1a1a1a1a1b1c2	ACT5041
																														O1a1a1a1a1a1b1d	BY135127
																														O1a1a1a1a1a1b1e	BY11120
																															O1a1a1a1a1a1b1e1	Y146786
																												O1a1a1a1a1a1c	Y31266
																													O1a1a1a1a1a1c1	Y31261
																														O1a1a1a1a1a1c1a	Y106084
																															O1a1a1a1a1a1c1a1	A23677
																													O1a1a1a1a1a1c2	Y142811
																														O1a1a1a1a1a1c2a	Y142799
																													O1a1a1a1a1a1c3	A22477
																												O1a1a1a1a1a1d	A12442
																													O1a1a1a1a1a1d1	A12441
																														O1a1a1a1a1a1d1a	Y76085
																															O1a1a1a1a1a1d1a1	F18517
																															O1a1a1a1a1a1d1a2	MF15781
																													O1a1a1a1a1a1d2	ACT1103
																												O1a1a1a1a1a1e	MF1071
																													O1a1a1a1a1a1e1	MF1074
																													O1a1a1a1a1a1e2	Y151570
																														O1a1a1a1a1a1e2a	Y151577
																															O1a1a1a1a1a1e2a1	MF37645
																												O1a1a1a1a1a1f	Y137973
																											O1a1a1a1a1a2	CTS4585
																												O1a1a1a1a1a2a	FGC15402
																														O1a1a1a1a1a2a1a	FGC15421
																													O1a1a1a1a1a2a2	A22480
																													O1a1a1a1a1a2a3	MF6395
																												O1a1a1a1a1a2b	MF6441
																										O1a1a1a1a1b	MF6080
																											O1a1a1a1a1b1	MF6083
																												O1a1a1a1a1b1a	MF6084
																									O1a1a1a1a2	MF1075
																										O1a1a1a1a2a	Y137919
																											O1a1a1a1a2a1	Y137924
																												O1a1a1a1a2a1a	MF21015
																												O1a1a1a1a2a1b	A23678
																									O1a1a1a1a3	ACT263
																										O1a1a1a1a3a	ACT133
																							O1a1a1a2	YP4610
																								O1a1a1a2a	AM00330
																									O1a1a1a2a1	AM00333
																										O1a1a1a2a1a	B388
																								O1a1a1a2b	SK1555
																							O1a1a1a3	L316.2
																								O1a1a1a3a	MF14293
																							O1a1a1a4	SK1530
																						O1a1a1b	F5498
																							O1a1a1b1	Z23406
																								O1a1a1b1a	M101
																									O1a1a1b1a1	SK1573
																									O1a1a1b1a2	A23875
																								O1a1a1b1b	Z23392
																									O1a1a1b1b1	Z23442
																										O1a1a1b1b1a	F22870
																											O1a1a1b1b1a1	SK1571
																												O1a1a1b1b1a1a	F20634
																										O1a1a1b1b1b	Z38638
																											O1a1a1b1b1b1	Z38636
																					O1a1a2	CTS8423
																						O1a1a2a	CTS52
																							O1a1a2a1	CTS701
																									O1a1a2a1a1	CTS716
																										O1a1a2a1a1a	MF2375
																										O1a1a2a1a1b	Y147549
																									O1a1a2a1a2	CTS352
																										O1a1a2a1a2a	CTS3144
																									O1a1a2a1a3	Z23271
																										O1a1a2a1a3a	Z7773
																										O1a1a2a1a3b	Z23326
																								O1a1a2a1b	Y157651
																									O1a1a2a1b1	CTS1992
																										O1a1a2a1b1a	CTS11040
																							O1a1a2a2	FGC66104
																								O1a1a2a2a	FGC66056
																									O1a1a2a2a1	FGC66089
																										O1a1a2a2a1a	FGC66085
																							O1a1a2a3	BY47757
																						O1a1a2b	F2444
																							O1a1a2b1	Y137055
																							O1a1a2b2	F1056
																								O1a1a2b2a	MF6458
																						O1a1a2c	SK1522
																							O1a1a2c1	MF6151
																				O1a1b	CTS5726
																					O1a1b1	SK1590
																						O1a1b1a	F1476
																							O1a1b1a1	MF2387
																								O1a1b1a1a	MF2958
																							O1a1b1a2	F2818
																						O1a1b1b	B390
																			O1a2	M110
																				O1a2a	F3288
																					O1a2a1	B392
																						O1a2a1a	B393
																					O1a2a2	F1600
																					O1a2a3	Y33185
																					O1a2a4	Z38625
																			O1a3	F1036
																				O1a3a	F1009
																					O1a3a1	F970
																		O1b	M268
																			O1b1	F2320
																				O1b1a	M1470
																					O1b1a1	PK4
																						O1b1a1a	M95
																							O1b1a1a1	F1803
																								O1b1a1a1a	F1252
																									O1b1a1a1a1	F2924
																										O1b1a1a1a1a	F2524
																											O1b1a1a1a1a1	M111
																												O1b1a1a1a1a1a	F2758
																													O1b1a1a1a1a1a1	Z24083
																														O1b1a1a1a1a1a1a	Z24089
																															O1b1a1a1a1a1a1a1	F923
																																O1b1a1a1a1a1a1a1a	A15721
																																	O1b1a1a1a1a1a1a1a1	CTS2022
																																			O1b1a1a1a1a1a1a1a1a1	F1399
																																				O1b1a1a1a1a1a1a1a1a1a	F2415
																																			O1b1a1a1a1a1a1a1a1a2	FGC79907
																																				O1b1a1a1a1a1a1a1a1a2a	FGC79915
																																	O1b1a1a1a1a1a1a1a2	Z24131
																															O1b1a1a1a1a1a1a2	SK1627
																																O1b1a1a1a1a1a1a2a	Z39410
																																O1b1a1a1a1a1a1a2b	Y26367
																														O1b1a1a1a1a1a1b	Z24088
																															O1b1a1a1a1a1a1b1	F19221
																																O1b1a1a1a1a1a1b1a	F19998
																												O1b1a1a1a1a1b	F2890
																													O1b1a1a1a1a1b1	Z24048
																														O1b1a1a1a1a1b1a	Z24050
																															O1b1a1a1a1a1b1a1	ACT523
																																O1b1a1a1a1a1b1a1a	F15446
																																	O1b1a1a1a1a1b1a1a1	F16217
																																O1b1a1a1a1a1b1a1b	ACT505
																																	O1b1a1a1a1a1b1a1b1	ACT553
																														O1b1a1a1a1a1b1b	F18990
																															O1b1a1a1a1a1b1b1	F4101
																													O1b1a1a1a1a1b2	Z24014
																										O1b1a1a1a1b	CTS5854
																											O1b1a1a1a1b1	Z23810
																												O1b1a1a1a1b1a	CTS7399
																													O1b1a1a1a1b1a1	FGC19713
																														O1b1a1a1a1b1a1a	FGC19707
																															O1b1a1a1a1b1a1a1	Z23849
																																O1b1a1a1a1b1a1a1a	Z23863
																																	O1b1a1a1a1b1a1a1a1	Y148916
																																		O1b1a1a1a1b1a1a1a1a	MF16225
																																			O1b1a1a1a1b1a1a1a1a1	Z46424
																																		O1b1a1a1a1b1a1a1a1b	BY73069
																																			O1b1a1a1a1b1a1a1a1b1	MF14663
																																				O1b1a1a1a1b1a1a1a1b1a	MF15901
																																		O1b1a1a1a1b1a1a1a1c	BY70435
																																			O1b1a1a1a1b1a1a1a1c1	MF16254
																																				O1b1a1a1a1b1a1a1a1c1a	F9934
																																	O1b1a1a1a1b1a1a1a2	F20080
																																		O1b1a1a1a1b1a1a1a2a	Z46431
																																		O1b1a1a1a1b1a1a1a2b	F15867
																																			O1b1a1a1a1b1a1a1a2b1	MF14552
																																	O1b1a1a1a1b1a1a1a3	FGC61038
																															O1b1a1a1a1b1a1a2	FGC19716
																																O1b1a1a1a1b1a1a2a	FGC14264
																																O1b1a1a1a1b1a1a2b	BY65986
																																	O1b1a1a1a1b1a1a2b1	MF14236
																												O1b1a1a1a1b1b	CTS651
																													O1b1a1a1a1b1b1	CTS9884
																											O1b1a1a1a1b2	Z23781
																												O1b1a1a1a1b2a	F4229
																													O1b1a1a1a1b2a1	F809
																														O1b1a1a1a1b2a1a	Z23795
																															O1b1a1a1a1b2a1a1	F2517
																																O1b1a1a1a1b2a1a1a	Z23790
																																O1b1a1a1a1b2a1a1b	MF2899
																															O1b1a1a1a1b2a1a2	A22937
																																O1b1a1a1a1b2a1a2a	A23595
																														O1b1a1a1a1b2a1b	FGC61076
																									O1b1a1a1a2	SK1630
																										O1b1a1a1a2a	ACT5539
																											O1b1a1a1a2a1	SK1636
																								O1b1a1a1b	F789
																									O1b1a1a1b1	FGC29900
																										O1b1a1a1b1a	B426
																											O1b1a1a1b1a1	FGC29907
																											O1b1a1a1b1a2	B427
																											O1b1a1a1b1a3	B424
																												O1b1a1a1b1a3a	B425
																											O1b1a1a1b1a4	F15214
																												O1b1a1a1b1a4a	F21137
																										O1b1a1a1b1b	Z39485
																											O1b1a1a1b1b1	F26580
																												O1b1a1a1b1b1a	F18676
																										O1b1a1a1b1c	B418
																											O1b1a1a1b1c1	B421
																											O1b1a1a1b1c2	F14906
																												O1b1a1a1b1c2a	F16164
																									O1b1a1a1b2	A22938
																										O1b1a1a1b2a	SK1646
																							O1b1a1a2	CTS350
																						O1b1a1b	F838
																							O1b1a1b1	CTS3857
																								O1b1a1b1a	F714
																									O1b1a1b1a1	F3357
																										O1b1a1b1a1a	F977
																											O1b1a1b1a1a1	F1199
																										O1b1a1b1a1b	CTS5664
																											O1b1a1b1a1b1	Y155291
																											O1b1a1b1a1b2	CTS2452
																								O1b1a1b1b	F15640
																									O1b1a1b1b1	F16212
																										O1b1a1b1b1a	F21414
																											O1b1a1b1b1a1	F15452
																					O1b1a2	Page59
																						O1b1a2a	F993
																							O1b1a2a1	F1759
																								O1b1a2a1a	CTS1127
																						O1b1a2b	F417
																							O1b1a2b1	F840
																								O1b1a2b1a	F1127
																							O1b1a2b2	CTS1451
																						O1b1a2c	CTS9996
																			O1b2	P49
																				O1b2a	F1942
																					O1b2a1	CTS9259
																						O1b2a1a	F1204
																							O1b2a1a1	CTS713
																								O1b2a1a1a	CTS1875
																									O1b2a1a1a1	CTS10682
																								O1b2a1a1b	Z24598
																								O1b2a1a1c	CTS203
																							O1b2a1a2	F2868
																								O1b2a1a2a	L682
																									O1b2a1a2a1	CTS723
																										O1b2a1a2a1a	CTS7620
																										O1b2a1a2a1b	A12446
																											O1b2a1a2a1b1	PH40
																								O1b2a1a2b	F940
																							O1b2a1a3	CTS10687
																								O1b2a1a3a	CTS1215
																						O1b2a1b	CTS562
																	O2	M122
																		O2a	M324
																			O2a1	L127.1
																				O2a1a	F1876
																					O2a1a1	F2159
																						O2a1a1a	F1867
																							O2a1a1a1	F852
																								O2a1a1a1a	F2266
																									O2a1a1a1a1	L599
																										O2a1a1a1a1a	Z24958
																											O2a1a1a1a1a1	Z43961
																												O2a1a1a1a1a1a	Z43963
																											O2a1a1a1a1a2	MF3005
																												O2a1a1a1a1a2a	MF2718
																											O2a1a1a1a1a3	MF14378
																										O2a1a1a1a1b	MF21009
																											O2a1a1a1a1b1	MF21989
																									O2a1a1a1a2	CTS3500
																								O2a1a1a1b	F854
																									O2a1a1a1b1	CTS1290
																										O2a1a1a1b1a	Z43966
																									O2a1a1a1b2	MF6584
																										O2a1a1a1b2a	F2687
																											O2a1a1a1b2a1	F2283
																							O2a1a1a2	MF14217
																								O2a1a1a2a	MF14214
																									O2a1a1a2a1	MF14184
																									O2a1a1a2a2	MF14337
																						O2a1a1b	F915
																							O2a1a1b1	F1478
																								O2a1a1b1a	PF5390
																									O2a1a1b1a1	CTS1936
																										O2a1a1b1a1a	Z43975
																										O2a1a1b1a1b	PH954
																											O2a1a1b1a1b1	MF14123
																									O2a1a1b1a2	FGC33994
																										O2a1a1b1a2a	A18028
																											O2a1a1b1a2a1	FGC34020
																												O2a1a1b1a2a1a	MF15799
																													O2a1a1b1a2a1a1	MF14471
																													O2a1a1b1a2a1a2	MF15777
																													O2a1a1b1a2a1a3	MF87544
																														O2a1a1b1a2a1a3a	MF14300
																										O2a1a1b1a2b	MF15143
																											O2a1a1b1a2b1	MF14383
																												O2a1a1b1a2b1a	MF16404
																											O2a1a1b1a2b2	ALK841.2
																												O2a1a1b1a2b2a	MF14702
																													O2a1a1b1a2b2a1	MF14599
																												O2a1a1b1a2b2b	MF14261
																													O2a1a1b1a2b2b1	MF14133
																													O2a1a1b1a2b2b2	MF14741
																														O2a1a1b1a2b2b2a	MF15713
																													O2a1a1b1a2b2b3	MF15080
																														O2a1a1b1a2b2b3a	MF14446
																													O2a1a1b1a2b2b4	MF15193
																														O2a1a1b1a2b2b4a	MF22423
																										O2a1a1b1a2c	MF14358
																									O2a1a1b1a3	MF2890
																								O2a1a1b1b	F905
																							O2a1a1b2	F4080
																								O2a1a1b2a	F4180
																								O2a1a1b2b	MF14286
																				O2a1b	IMS-JST002611
																					O2a1b1	F18
																						O2a1b1a	FGC12511
																							O2a1b1a1	F117
																								O2a1b1a1a	CTS4598
																									O2a1b1a1a1	CTS4263
																										O2a1b1a1a1a	F11
																											O2a1b1a1a1a1	F325
																												O2a1b1a1a1a1a	F632
																													O2a1b1a1a1a1a1	F110
																														O2a1b1a1a1a1a1a	CTS3663
																															O2a1b1a1a1a1a1a1	F17
																																O2a1b1a1a1a1a1a1a	F377
																																	O2a1b1a1a1a1a1a1a1	CTS7789
																																		O2a1b1a1a1a1a1a1a1a	F856
																																			O2a1b1a1a1a1a1a1a1a1	F1418
																																			O2a1b1a1a1a1a1a1a1a2	Z25098
																																				O2a1b1a1a1a1a1a1a1a2a	Z25097
																																					O2a1b1a1a1a1a1a1a1a2a1	Z7776.2
																																						O2a1b1a1a1a1a1a1a1a2a1a	MF16427
																																					O2a1b1a1a1a1a1a1a1a2a2	Z46463
																																		O2a1b1a1a1a1a1a1a1b	MF14538
																																			O2a1b1a1a1a1a1a1a1b1	MF15136
																																	O2a1b1a1a1a1a1a1a2	CTS7501
																																		O2a1b1a1a1a1a1a1a2a	CTS1621
																																		O2a1b1a1a1a1a1a1a2b	BY57051
																																O2a1b1a1a1a1a1a1b	Z43876
																																	O2a1b1a1a1a1a1a1b1	F793
																																		O2a1b1a1a1a1a1a1b1a	Y141213
																																			O2a1b1a1a1a1a1a1b1a1	MF7453
																																		O2a1b1a1a1a1a1a1b1b	Y65006
																																			O2a1b1a1a1a1a1a1b1b1	A22940
																																		O2a1b1a1a1a1a1a1b1c	MF30934
																																			O2a1b1a1a1a1a1a1b1c1	MF31684
																																		O2a1b1a1a1a1a1a1b1d	Y150436
																																		O2a1b1a1a1a1a1a1b1e	CTS10514
																																		O2a1b1a1a1a1a1a1b1f	F931
																														O2a1b1a1a1a1a1b	Z43879
																															O2a1b1a1a1a1a1b1	Y20951
																																O2a1b1a1a1a1a1b1a	Y20939
																																	O2a1b1a1a1a1a1b1a1	Y20932
																																		O2a1b1a1a1a1a1b1a1a	FGC66016
																																	O2a1b1a1a1a1a1b1a2	MF6724
																																O2a1b1a1a1a1a1b1b	FGC23704
																															O2a1b1a1a1a1a1b2	CTS7271
																														O2a1b1a1a1a1a1c	MF6739
																															O2a1b1a1a1a1a1c1	MF6754
																																O2a1b1a1a1a1a1c1a	FGC1419
																													O2a1b1a1a1a1a2	F16340
																														O2a1b1a1a1a1a2a	F15504
																												O2a1b1a1a1a1b	F38
																													O2a1b1a1a1a1b1	F386
																														O2a1b1a1a1a1b1a	A22321
																															O2a1b1a1a1a1b1a1	F136
																																O2a1b1a1a1a1b1a1a	F2405
																																	O2a1b1a1a1a1b1a1a1	A22322
																																O2a1b1a1a1a1b1a1b	MF6807
																																	O2a1b1a1a1a1b1a1b1	FGC33600
																																		O2a1b1a1a1a1b1a1b1a	MF6812
																																			O2a1b1a1a1a1b1a1b1a1	MF6813
																																				O2a1b1a1a1a1b1a1b1a1a	A7966
																																					O2a1b1a1a1a1b1a1b1a1a1	MF15642
																																						O2a1b1a1a1a1b1a1b1a1a1a	BY42738
																																				O2a1b1a1a1a1b1a1b1a1b	MF7684
																																			O2a1b1a1a1a1b1a1b1a2	MF7764
																																				O2a1b1a1a1a1b1a1b1a2a	MF15180
																																	O2a1b1a1a1a1b1a1b2	MF14988
																																		O2a1b1a1a1a1b1a1b2a	MF46495
																																	O2a1b1a1a1a1b1a1b3	F2091
																																O2a1b1a1a1a1b1a1c	MF7704
																																	O2a1b1a1a1a1b1a1c1	Z46464
																																		O2a1b1a1a1a1b1a1c1a	MF14553
																																		O2a1b1a1a1a1b1a1c1b	MF14959
																																	O2a1b1a1a1a1b1a1c2	MF7673
																																		O2a1b1a1a1a1b1a1c2a	MF7721
																																O2a1b1a1a1a1b1a1d	MF7678
																																	O2a1b1a1a1a1b1a1d1	MF7716
																																		O2a1b1a1a1a1b1a1d1a	BY68984
																																O2a1b1a1a1a1b1a1e	MF7664
																																O2a1b1a1a1a1b1a1f	A22623
																															O2a1b1a1a1a1b1a2	MF6816
																																O2a1b1a1a1a1b1a2a	MF6826
																																	O2a1b1a1a1a1b1a2a1	MF6836
																															O2a1b1a1a1a1b1a3	MF14171
																														O2a1b1a1a1a1b1b	A10289
																													O2a1b1a1a1a1b2	MF7683
																												O2a1b1a1a1a1c	F12
																													O2a1b1a1a1a1c1	F319
																														O2a1b1a1a1a1c1a	F683
																															O2a1b1a1a1a1c1a1	Z25153
																														O2a1b1a1a1a1c1b	MF6882
																															O2a1b1a1a1a1c1b1	MF6883
																																O2a1b1a1a1a1c1b1a	MF7876
																													O2a1b1a1a1a1c2	F16631
																												O2a1b1a1a1a1d	F930
																													O2a1b1a1a1a1d1	F2685
																														O2a1b1a1a1a1d1a	F2202
																												O2a1b1a1a1a1e	F1365
																													O2a1b1a1a1a1e1	Y15976
																														O2a1b1a1a1a1e1a	Y16154
																															O2a1b1a1a1a1e1a1	Y26383
																																O2a1b1a1a1a1e1a1a	Y26386
																																	O2a1b1a1a1a1e1a1a1	SK1686
																																		O2a1b1a1a1a1e1a1a1a	Y46258
																																		O2a1b1a1a1a1e1a1a1b	F15966
																																O2a1b1a1a1a1e1a1b	MF6915
																																O2a1b1a1a1a1e1a1c	MF8048
																															O2a1b1a1a1a1e1a2	MF6932
																																O2a1b1a1a1a1e1a2a	MF6934
																																	O2a1b1a1a1a1e1a2a1	MF6935
																																	O2a1b1a1a1a1e1a2a2	MF6938
																															O2a1b1a1a1a1e1a3	MF6974
																															O2a1b1a1a1a1e1a4	MF6953
																														O2a1b1a1a1a1e1b	F1200
																													O2a1b1a1a1a1e2	FGC54486
																														O2a1b1a1a1a1e2a	MF2463
																															O2a1b1a1a1a1e2a1	FGC54507
																																O2a1b1a1a1a1e2a1a	FGC54505
																																	O2a1b1a1a1a1e2a1a1	MF8087
																															O2a1b1a1a1a1e2a2	MF2473
																																O2a1b1a1a1a1e2a2a	ACT2777
																																	O2a1b1a1a1a1e2a2a1	ACT2793
																																O2a1b1a1a1a1e2a2b	MF2464
																														O2a1b1a1a1a1e2b	MF1083
																															O2a1b1a1a1a1e2b1	MF9483
																																O2a1b1a1a1a1e2b1a	MF9486
																														O2a1b1a1a1a1e2c	Y163866
																														O2a1b1a1a1a1e2d	MF8460
																												O2a1b1a1a1a1f	CTS12877
																													O2a1b1a1a1a1f1	F2527
																														O2a1b1a1a1a1f1a	CTS5409
																															O2a1b1a1a1a1f1a1	CTS9711
																																O2a1b1a1a1a1f1a1a	MF1309
																																O2a1b1a1a1a1f1a1b	CTS238
																															O2a1b1a1a1a1f1a2	MF1085
																																O2a1b1a1a1a1f1a2a	F16712
																																	O2a1b1a1a1a1f1a2a1	MF8589
																																		O2a1b1a1a1a1f1a2a1a	MF8504
																														O2a1b1a1a1a1f1b	F2941
																															O2a1b1a1a1a1f1b1	MF21873
																													O2a1b1a1a1a1f2	CTS2209
																														O2a1b1a1a1a1f2a	MF17282
																														O2a1b1a1a1a1f2b	CTS4241
																												O2a1b1a1a1a1g	F723
																													O2a1b1a1a1a1g1	ACT5910
																												O2a1b1a1a1a1h	CTS2154
																													O2a1b1a1a1a1h1	CTS2107
																														O2a1b1a1a1a1h1a	MF8526
																															O2a1b1a1a1a1h1a1	MF8582
																															O2a1b1a1a1a1h1a2	Y135777
																														O2a1b1a1a1a1h1b	Y150148
																													O2a1b1a1a1a1h2	MF1094
																														O2a1b1a1a1a1h2a	MF8512
																												O2a1b1a1a1a1i	SK1691
																													O2a1b1a1a1a1ia	SK1692
																													O2a1b1a1a1a1ib	MF20736
																														O2a1b1a1a1a1ib1	MF20650
																												O2a1b1a1a1a1j	MF1107
																													O2a1b1a1a1a1j1	MF1145
																														O2a1b1a1a1a1j1a	MF8613
																												O2a1b1a1a1a1k	SK1676
																												O2a1b1a1a1a1l	FGC60773
																													O2a1b1a1a1a1l1	FGC60796
																											O2a1b1a1a1a2	MF6982
																										O2a1b1a1a1b	PH203
																											O2a1b1a1a1b1	FGC1803.2
																												O2a1b1a1a1b1a	MF7030
																												O2a1b1a1a1b1b	MF54301
																													O2a1b1a1a1b1b1	MF55477
																											O2a1b1a1a1b2	MF7062
																									O2a1b1a1a2	MF7120
																								O2a1b1a1b	MF1237
																									O2a1b1a1b1	MF1292
																							O2a1b1a2	F449
																								O2a1b1a2a	F238
																									O2a1b1a2a1	F271
																										O2a1b1a2a1a	F134
																											O2a1b1a2a1a1	CTS10846
																												O2a1b1a2a1a1a	CTS6279
																													O2a1b1a2a1a1a1	F25545
																														O2a1b1a2a1a1a1a	MF1150
																														O2a1b1a2a1a1a1b	CTS9595
																															O2a1b1a2a1a1a1b1	F20523
																														O2a1b1a2a1a1a1c	MF8759
																														O2a1b1a2a1a1a1d	MF18144
																															O2a1b1a2a1a1a1d1	SK1698
																																O2a1b1a2a1a1a1d1a	F18567
																														O2a1b1a2a1a1a1e	MF17397
																													O2a1b1a2a1a1a2	MF14479
																														O2a1b1a2a1a1a2a	F1435
																														O2a1b1a2a1a1a2b	MF14140
																												O2a1b1a2a1a1b	MF7188
																													O2a1b1a2a1a1b1	MF7207
																											O2a1b1a2a1a2	F724
																												O2a1b1a2a1a2a	F288
																												O2a1b1a2a1a2b	MF8764
																											O2a1b1a2a1a3	F21824
																												O2a1b1a2a1a3a	F18913
																													O2a1b1a2a1a3a1	MF7243
																														O2a1b1a2a1a3a1a	MF8757
																														O2a1b1a2a1a3a1b	MF7251
																															O2a1b1a2a1a3a1b1	MF7250
																																O2a1b1a2a1a3a1b1a	MF7254
																																	O2a1b1a2a1a3a1b1a1	MF15031
																																	O2a1b1a2a1a3a1b1a2	MF7248
																																		O2a1b1a2a1a3a1b1a2a	MF7285
																												O2a1b1a2a1a3b	MF1156
																													O2a1b1a2a1a3b1	MF1159
																													O2a1b1a2a1a3b2	BY38244
																														O2a1b1a2a1a3b2a	MF8936
																														O2a1b1a2a1a3b2b	MF14440
																												O2a1b1a2a1a3c	MF2694
																											O2a1b1a2a1a4	MF14276
																												O2a1b1a2a1a4a	MF14692
																									O2a1b1a2a2	CTS6840
																										O2a1b1a2a2a	MF14314
																											O2a1b1a2a2a1	MF17568
																												O2a1b1a2a2a1a	MF14621
																												O2a1b1a2a2a1b	MF15708
																													O2a1b1a2a2a1b1	MF16703
																														O2a1b1a2a2a1b1a	MF16618
																										O2a1b1a2a2b	CTS679
																											O2a1b1a2a2b1	CTS5105
																									O2a1b1a2a3	MF17449
																								O2a1b1a2b	F1266
																									O2a1b1a2b1	MF14122
																									O2a1b1a2b2	F1403
																										O2a1b1a2b2a	MF9185
																											O2a1b1a2b2a1	MF9226
																										O2a1b1a2b2b	F2459
																											O2a1b1a2b2b1	F2725
																												O2a1b1a2b2b1a	MF7359
																													O2a1b1a2b2b1a1	MF14177
																														O2a1b1a2b2b1a1a	MF9266
																													O2a1b1a2b2b1a2	MF9267
																												O2a1b1a2b2b1b	MF7360
																													O2a1b1a2b2b1b1	BY139081
																													O2a1b1a2b2b1b2	MF14405
																						O2a1b1b	CTS498
																							O2a1b1b1	MF1348
																							O2a1b1b2	CTS5907
																					O2a1b2	FGC3750
																						O2a1b2a	CTS879
																							O2a1b2a1	SK1670
																								O2a1b2a1a	SK1671
																									O2a1b2a1a1	M164
																								O2a1b2a1b	MF9326
																							O2a1b2a2	Z38921
																								O2a1b2a2a	Z38920
																			O2a2	IMS-JST021354
																				O2a2a	M188
																					O2a2a1	F2588
																						O2a2a1a	CTS445
																							O2a2a1a1	CTS201
																								O2a2a1a1a	M159
																									O2a2a1a1a1	CTS1763
																										O2a2a1a1a1a	CTS1980
																										O2a2a1a1a1b	Z35182.2
																											O2a2a1a1a1b1	MF2541
																												O2a2a1a1a1b1a	MF2545
																													O2a2a1a1a1b1a1	MF2546
																										O2a2a1a1a1c	MF2558
																											O2a2a1a1a1c1	MF2549
																												O2a2a1a1a1c1a	MF2548
																								O2a2a1a1b	FGC50625
																									O2a2a1a1b1	Y52580.2
																										O2a2a1a1b1a	A11476.2
																									O2a2a1a1b2	FGC50672
																										O2a2a1a1b2a	Y172144
																										O2a2a1a1b2b	FGC50710
																											O2a2a1a1b2b1	PH2141
																											O2a2a1a1b2b2	FGC50673
																												O2a2a1a1b2b2a	MF15268
																													O2a2a1a1b2b2a1	FGC32772
																														O2a2a1a1b2b2a1a	MF18626
																														O2a2a1a1b2b2a1b	MF15693
																															O2a2a1a1b2b2a1b1	MF15348
																												O2a2a1a1b2b2b	SK1702
																													O2a2a1a1b2b2b1	FGC50535
																														O2a2a1a1b2b2b1a	FGC50540
																															O2a2a1a1b2b2b1a1	FGC50613
																																O2a2a1a1b2b2b1a1a	FGC50638
																															O2a2a1a1b2b2b1a2	BY46774
																															O2a2a1a1b2b2b1a3	Z46466
																																O2a2a1a1b2b2b1a3a	Z46467
																															O2a2a1a1b2b2b1a4	MF15589
																																O2a2a1a1b2b2b1a4a	Z46468
																															O2a2a1a1b2b2b1a5	MF3726
																														O2a2a1a1b2b2b1b	Y65106.2
																													O2a2a1a1b2b2b2	BY54719
																													O2a2a1a1b2b2b3	Y165598
																													O2a2a1a1b2b2b4	BY63584
																														O2a2a1a1b2b2b4a	BY55062
																							O2a2a1a2	M7
																								O2a2a1a2a	CTS10944
																									O2a2a1a2a1	F1276
																										O2a2a1a2a1a	CTS6489
																											O2a2a1a2a1a1	F1275
																												O2a2a1a2a1a1a	Z25411
																													O2a2a1a2a1a1a1	F827
																														O2a2a1a2a1a1a1a	M113
																													O2a2a1a2a1a1a2	Z25398
																														O2a2a1a2a1a1a2a	F1100
																															O2a2a1a2a1a1a2a1	F1234
																																O2a2a1a2a1a1a2a1a	F1411
																																	O2a2a1a2a1a1a2a1a1	N5
																																		O2a2a1a2a1a1a2a1a1a	F15848
																																			O2a2a1a2a1a1a2a1a1a1	F22246
																																				O2a2a1a2a1a1a2a1a1a1a	F16501
																																O2a2a1a2a1a1a2a1b	MF3299
																															O2a2a1a2a1a1a2a2	Y161282
																														O2a2a1a2a1a1a2b	Z25400
																												O2a2a1a2a1a1b	CTS6579
																													O2a2a1a2a1a1b1	F14832
																														O2a2a1a2a1a1b1a	F15414
																													O2a2a1a2a1a1b2	CTS123
																												O2a2a1a2a1a1c	F14807
																													O2a2a1a2a1a1c1	F15259
																														O2a2a1a2a1a1c1a	F19315
																										O2a2a1a2a1b	F1863
																											O2a2a1a2a1b1	F1134
																												O2a2a1a2a1b1a	F1262
																									O2a2a1a2a2	Y26403
																									O2a2a1a2a3	MF9858
																										O2a2a1a2a3a	MF21594
																						O2a2a1b	F1837
																					O2a2a2	F862
																						O2a2a2a	F879
																							O2a2a2a1	F1226
																								O2a2a2a1a	F2937
																									O2a2a2a1a1	F4125
																										O2a2a2a1a1a	F2859
																											O2a2a2a1a1a1	F2696
																							O2a2a2a2	CTS8490
																				O2a2b	P164
																					O2a2b1	M134
																						O2a2b1a	F450
																							O2a2b1a1	M117
																								O2a2b1a1a	M133
																									O2a2b1a1a1	F8
																										O2a2b1a1a1a	A9459
																											O2a2b1a1a1a1	F438
																												O2a2b1a1a1a1a	F14366
																													O2a2b1a1a1a1a1	Y17728
																														O2a2b1a1a1a1a1a	F316
																															O2a2b1a1a1a1a1a1	F155
																																	O2a2b1a1a1a1a1a1a1	F813
																																		O2a2b1a1a1a1a1a1a1a	Y20928
																																			O2a2b1a1a1a1a1a1a1a1	MF11011
																																			O2a2b1a1a1a1a1a1a1a2	Y137613
																																				O2a2b1a1a1a1a1a1a1a2a	Y137637
																																			O2a2b1a1a1a1a1a1a1a3	MF21
																																			O2a2b1a1a1a1a1a1a1a4	MF2610
																																			O2a2b1a1a1a1a1a1a1a5	CTS7181
																																		O2a2b1a1a1a1a1a1a1b	Y154488
																																			O2a2b1a1a1a1a1a1a1b1	MF1
																																				O2a2b1a1a1a1a1a1a1b1a	FGC29526
																																					O2a2b1a1a1a1a1a1a1b1a1	MF3
																																			O2a2b1a1a1a1a1a1a1b2	MF4344
																																		O2a2b1a1a1a1a1a1a1c	Y138426
																																		O2a2b1a1a1a1a1a1a1d	MF4
																																		O2a2b1a1a1a1a1a1a1e	MF5
																																		O2a2b1a1a1a1a1a1a1f	MF6
																																		O2a2b1a1a1a1a1a1a1g	F22162
																																O2a2b1a1a1a1a1a1b	Y138597
																																	O2a2b1a1a1a1a1a1b1	Y146618
																																O2a2b1a1a1a1a1a1c	MF3302
																																O2a2b1a1a1a1a1a1d	F1494
																															O2a2b1a1a1a1a1a2	Y138400
																														O2a2b1a1a1a1a1b	ACT7481
																																O2a2b1a1a1a1a1b1a	F2137
																																		O2a2b1a1a1a1a1b1a1a	F1442
																																			O2a2b1a1a1a1a1b1a1a1	F2265
																																				O2a2b1a1a1a1a1b1a1a1a	F1123
																																					O2a2b1a1a1a1a1b1a1a1a1	F1369
																																			O2a2b1a1a1a1a1b1a1a2	BY79235
																																		O2a2b1a1a1a1a1b1a1b	MF15397
																																			O2a2b1a1a1a1a1b1a1b1	A9511.2
																																				O2a2b1a1a1a1a1b1a1b1a	MF14333
																																					O2a2b1a1a1a1a1b1a1b1a1	MF14879
																																			O2a2b1a1a1a1a1b1a1b2	MF15764
																																				O2a2b1a1a1a1a1b1a1b2a	MF15281
																																					O2a2b1a1a1a1a1b1a1b2a1	MF14573
																																					O2a2b1a1a1a1a1b1a1b2a2	MF15326
																																				O2a2b1a1a1a1a1b1a1b2b	MF14287
																																	O2a2b1a1a1a1a1b1a2	A16636
																																		O2a2b1a1a1a1a1b1a2a	A16639
																																			O2a2b1a1a1a1a1b1a2a1	BY78763
																																				O2a2b1a1a1a1a1b1a2a1a	MF14664
																																		O2a2b1a1a1a1a1b1a2b	MF14296
																																	O2a2b1a1a1a1a1b1a3	BY71124
																															O2a2b1a1a1a1a1b2	Y170157
																																O2a2b1a1a1a1a1b2a	Z46469
																																	O2a2b1a1a1a1a1b2a1	Z46473
																																		O2a2b1a1a1a1a1b2a1a	Z46475
																														O2a2b1a1a1a1a1c	Z25907
																															O2a2b1a1a1a1a1c1	Z25908
																															O2a2b1a1a1a1a1c2	MF14162
																													O2a2b1a1a1a1a2	F1074
																												O2a2b1a1a1a1b	MF3168
																											O2a2b1a1a1a2	FGC23469
																												O2a2b1a1a1a2a	FGC23471
																													O2a2b1a1a1a2a1	F666
																														O2a2b1a1a1a2a1a	F310
																															O2a2b1a1a1a2a1a1	F402
																																O2a2b1a1a1a2a1a1a	F1531
																											O2a2b1a1a1a3	F14249
																												O2a2b1a1a1a3a	Z25853
																													O2a2b1a1a1a3a1	CTS4789
																														O2a2b1a1a1a3a1a	CTS5492
																															O2a2b1a1a1a3a1a1	CTS6987
																																O2a2b1a1a1a3a1a1a	Z25902
																																	O2a2b1a1a1a3a1a1a1	Z42620
																															O2a2b1a1a1a3a1a2	F20963
																											O2a2b1a1a1a4	F14494
																												O2a2b1a1a1a4a	CTS4658
																													O2a2b1a1a1a4a1	CTS5308
																													O2a2b1a1a1a4a2	Z25928
																														O2a2b1a1a1a4a2a	Z25930
																															O2a2b1a1a1a4a2a1	SK1730
																																O2a2b1a1a1a4a2a1a	Z26030
																																	O2a2b1a1a1a4a2a1a1	Z26031
																																O2a2b1a1a1a4a2a1b	Z26010
																															O2a2b1a1a1a4a2a2	A9462
																															O2a2b1a1a1a4a2a3	Z39704
																																O2a2b1a1a1a4a2a3a	B456
																																	O2a2b1a1a1a4a2a3a1	Z39706
																											O2a2b1a1a1a5	F995
																												O2a2b1a1a1a5a	F2422
																										O2a2b1a1a1b	CTS7634
																											O2a2b1a1a1b1	F574
																												O2a2b1a1a1b1a	F474
																													O2a2b1a1a1b1a1	F317
																														O2a2b1a1a1b1a1a	F3039
																														O2a2b1a1a1b1a1b	Y29861
																											O2a2b1a1a1b2	CTS5488
																										O2a2b1a1a1c	CTS10738
																											O2a2b1a1a1c1	CTS7316
																												O2a2b1a1a1c1a	M1732
																													O2a2b1a1a1c1a1	CTS9678
																														O2a2b1a1a1c1a1a	CTS8848
																															O2a2b1a1a1c1a1a1	Z39663
																															O2a2b1a1a1c1a1a2	M1513
																																O2a2b1a1a1c1a1a2a	M1557
																												O2a2b1a1a1c1b	A9457
																													O2a2b1a1a1c1b1	F17158
																										O2a2b1a1a1d	YP4864
																											O2a2b1a1a1d1	Z44068
																												O2a2b1a1a1d1a	F5524
																													O2a2b1a1a1d1a1	F5525
																											O2a2b1a1a1d2	Z44073
																												O2a2b1a1a1d2a	Z44071
																										O2a2b1a1a1e	Z44091
																											O2a2b1a1a1e1	Z44092
																											O2a2b1a1a1e2	SK1740
																												O2a2b1a1a1e2a	Y137505
																													O2a2b1a1a1e2a1	Y137559
																										O2a2b1a1a1f	A7011.2
																										O2a2b1a1a1g	MF2689
																										O2a2b1a1a1h	Y138333
																								O2a2b1a1b	CTS4960
																							O2a2b1a2	F122
																								O2a2b1a2a	F114
																									O2a2b1a2a1	F79
																										O2a2b1a2a1a	F46
																											O2a2b1a2a1a1	FGC16847
																												O2a2b1a2a1a1a	F14184
																													O2a2b1a2a1a1a1	F48
																														O2a2b1a2a1a1a1a	CTS1011
																															O2a2b1a2a1a1a1a1	F55
																																O2a2b1a2a1a1a1a1a	F152
																																	O2a2b1a2a1a1a1a1a1	F700
																																		O2a2b1a2a1a1a1a1a1a	F14475
																																			O2a2b1a2a1a1a1a1a1a1	F15657
																																				O2a2b1a2a1a1a1a1a1a1a	F2505
																															O2a2b1a2a1a1a1a2	CTS3149
																													O2a2b1a2a1a1a2	F242
																														O2a2b1a2a1a1a2a	CTS6493
																															O2a2b1a2a1a1a2a1	F273
																																O2a2b1a2a1a1a2a1a	CTS7377
																																	O2a2b1a2a1a1a2a1a1	CTS4266
																																		O2a2b1a2a1a1a2a1a1a	Z26108
																																			O2a2b1a2a1a1a2a1a1a1	Y30135
																																				O2a2b1a2a1a1a2a1a1a1a	F2173
																																					O2a2b1a2a1a1a2a1a1a1a1	F1458
																												O2a2b1a2a1a1b	F2887
																													O2a2b1a2a1a1b1	CTS1346
																														O2a2b1a2a1a1b1a	F3607
																															O2a2b1a2a1a1b1a1	F3525
																														O2a2b1a2a1a1b1b	CTS3763
																															O2a2b1a2a1a1b1b1	A9472
																															O2a2b1a2a1a1b1b2	FGC16863
																																O2a2b1a2a1a1b1b2a	L1360
																																	O2a2b1a2a1a1b1b2a1	FGC16889
																																O2a2b1a2a1a1b1b2b	SK1768
																																	O2a2b1a2a1a1b1b2b1	F4249
																																		O2a2b1a2a1a1b1b2b1a	FGC23868
																																	O2a2b1a2a1a1b1b2b2	CTS2882
																																		O2a2b1a2a1a1b1b2b2a	CTS335
																											O2a2b1a2a1a2	CTS11910
																												O2a2b1a2a1a2a	CTS53
																													O2a2b1a2a1a2a1	CTS6373
																														O2a2b1a2a1a2a1a	A9473
																											O2a2b1a2a1a3	F3386
																											O2a2b1a2a1a4	Y29828
																												O2a2b1a2a1a4a	F735
																													O2a2b1a2a1a4a1	L932.2
																														O2a2b1a2a1a4a1a	FGC34973
																													O2a2b1a2a1a4a2	F1739
																									O2a2b1a2a2	F743
																										O2a2b1a2a2a	CTS8481.2
																											O2a2b1a2a2a1	CTS4325
																												O2a2b1a2a2a1a	A16629
																												O2a2b1a2a2a1b	CTS682
																										O2a2b1a2a2b	F748
																											O2a2b1a2a2b1	F1448
																												O2a2b1a2a2b1a	F728
																					O2a2b2	AM01822
																						O2a2b2a	AM01856
																							O2a2b2a1	N7
																								O2a2b2a1a	F4110
																									O2a2b2a1a1	F4068
																									O2a2b2a1a2	SK1780
																								O2a2b2a1b	F4124
																									O2a2b2a1b1	IMS-JST008425p6
																									O2a2b2a1b2	BY15188
																										O2a2b2a1b2a	F16411
																							O2a2b2a2	AM01845
																								O2a2b2a2a	F717
																									O2a2b2a2a1	F3612
																									O2a2b2a2a2	SK1783
																								O2a2b2a2b	AM01847
																									O2a2b2a2b1	A17418
																									O2a2b2a2b2	AM01756
																										O2a2b2a2b2a	B450
																										O2a2b2a2b2b	AM00472
																											O2a2b2a2b2b1	F18942
																										O2a2b2a2b2c	A16427
																						O2a2b2b	A16433
																							O2a2b2b1	A16438
																								O2a2b2b1a	SK1775
																									O2a2b2b1a1	SK1774
																								O2a2b2b1b	A16440
																		O2b	F742
																			O2b1	F1150
																				O2b1a	F837
																					O2b1a1	CTS4452
																						O2b1a1a	Y173835
																							O2b1a1a1	CTS6801
																								O2b1a1a1a	CTS798
																									O2b1a1a1a1	F1025
																								O2b1a1a1b	MF14654
																									O2b1a1a1b1	BY180079
																							O2b1a1a2	FGC62833
																					O2b1a2	BY37723
																				O2b1b	BY56049
																			O2b2	F1055
																				O2b2a	F3021
														K2	M526
															K2a	M2308
																M	SK1828
																	M1	M4
																		M1a	Z30983
																			M1a1	SK1846
																				M1a1a	BY8413.1
																					M1a1a2	Z30928
																						M1a1a2a	P87
																							M1a1a2a1	M104_1
																								M1a1a2a2a	M16
																								M1a1a2a2b	M83
																						M1a1a3b	B256
																							M1a1a3b1	Z33322
																								M1a1a3b1b	Z42297
																								M1a1a3b1c	B269
																							M1a1a3b2	Y26293
																	M2	M353
																		M2a	M177
																	M3	P117
																S	B254
																	S1	B255
																		S1a	Z41335
																			S1a1	Z42413
																				S1a1a	Z41513
																					S1a1a1	P308
																				S1a1b	M230
																					S1a1b1	M254
																						S1a1b1a	P57
																						S1a1b1b	P61
																						S1a1b1c	P83.1
																						S1a1b1d	SK1891
																							S1a1b1d1	Z41342
																								S1a1b1d1a	M226.1
																			S1a2	Z41767
																				S1a2a	P307
																		S1d	SK1806
																	S2	P378
															K2b	P331
																P or K2b2	P295
																	P1 or K2b2a	CTS196
																		Q or K2b2a1	M242
																		R or K2b2a2	M207
																	P2	F20148
																		P2b	F24883
																			P2b1	BY49600
																	Q	M242
																		Q1	F903
																			Q1a	F1096
																				Q1a1	F746
																					Q1a1a	M120
																						Q1a1a1	F1626
																							Q1a1a1a	F4528
																						Q1a1a2	M7417
																					Q1a1b	B143
																						Q1a1b2	B280
																							Q1a1b2b	B281
																				Q1a2	M25
																					Q1a2a	L712
																						Q1a2a1	L715
																			Q1b	M346
																				Q1b1	L53
																					Q1b1a	L54
																						Q1b1a1	CTS11969
																							Q1b1a1a	M3
																								Q1b1a1a1	CTS2610
																									Q1b1a1a1e	CTS11357
																										Q1b1a1a1e1	CTS11330
																											Q1b1a1a1e1a	CTS10359
																												Q1b1a1a1e1a1	CTS479
																										Q1b1a1a1e2	Z7796
																									Q1b1a1a1f	Z35840
																									Q1b1a1a1h	Z5906
																										Q1b1a1a1h1	CTS4000
																											Q1b1a1a1h1a	Z5907
																												Q1b1a1a1h1a4	B37
																													Q1b1a1a1h1a4b	B38
																									Q1b1a1a1i	Z5908
																										Q1b1a1a1i1	Z19318
																											Q1b1a1a1i1a	Z5910
																												Q1b1a1a1i1a1	Z5911
																													Q1b1a1a1i1a1a	Z5912
																												Q1b1a1a1i1a2	Z35921
																									Q1b1a1a1j	Y789
																										Q1b1a1a1j1	Z5914
																											Q1b1a1a1j1a	Y818
																									Q1b1a1a1k	Z19432
																										Q1b1a1a1k1	Z6658
																											Q1b1a1a1k1a	Z5916
																									Q1b1a1a1l	SK281
																									Q1b1a1a1m	CTS2731
																									Q1b1a1a1p	B43
																								Q1b1a1a2	FGC8469
																							Q1b1a1b	E324
																						Q1b1a2	CTS1780
																								Q1b1a2a1	Y2816
																										Q1b1a2a1a1	L569
																									Q1b1a2a1b	Z782
																								Q1b1a2a2	YP910
																						Q1b1a3	L330
																							Q1b1a3b	B287
																								Q1b1a3b1	BZ99
																									Q1b1a3b1a	B30
																					Q1b2a	F4674
																		Q2	L275
																			Q2a	F1213
																				Q2a1	L214
																					Q2a1a	L245
																						Q2a1a1	Y2202
																							Q2a1a1a	FGC1917
																								Q2a1a1a1	FGC1897
																									Q2a1a1a1a	FGC1929
																										Q2a1a1a1a1	FGC1933
																											Q2a1a1a1a1a	FGC1932
																												Q2a1a1a1a1a1	FGC4846
																												Q2a1a1a1a1a2	FGC1905
																													Q2a1a1a1a1a2a	FGC1898
																												Q2a1a1a1a1a3	YP1004
																													Q2a1a1a1a1a3a	YP3926
																													Q2a1a1a1a1a3b	BZ40
																										Q2a1a1a1a2	FGC23583
																											Q2a1a1a1a2a	FGC23580
																												Q2a1a1a1a2a1	YP1010
																													Q2a1a1a1a2a1a	YP1009
																												Q2a1a1a1a2a2	FGC23581
																														Q2a1a1a1a2a2a1	BZ53
																								Q2a1a1a2	YP730
																									Q2a1a1a2a	ACT296
																									Q2a1a1a2b	YP740
																						Q2a1a2	YP745
																							Q2a1a2a	YP1095
																								Q2a1a2a1	YP1096
																									Q2a1a2a1a	YP1240
																										Q2a1a2a1a1	YP1236
																											Q2a1a2a1a1a	YP1237
																												Q2a1a2a1a1a1	BZ3048
																													Q2a1a2a1a1a1a	BZ3043
																										Q2a1a2a1a2	SK433
																									Q2a1a2a1b	BZ6
																							Q2a1a2b	Y11199
																								Q2a1a2b1	FGC42712
																								Q2a1a2b2	Y11198
																									Q2a1a2b2a	FGC34074
																							Q2a1a2c	BZ278
																								Q2a1a2c1	BZ281
																									Q2a1a2c1a	BZ10
																										Q2a1a2c1a1	BZ3074
																									Q2a1a2c1b	BZ5047
																						Q2a1a3	BZ310
																							Q2a1a3a	BZ306
																								Q2a1a3a1	BZ3057
																									Q2a1a3a1b	Y166175
																						Q2a1a4	BZ3896
																							Q2a1a4a	BZ3900
																								Q2a1a4a1	L619.2
																								Q2a1a4a2	BZ3920
																					Q2a1b	BZ384
																					Q2a1c	FGC4607
																						Q2a1c1	FGC4626
																							Q2a1c1a	FGC4888
																							Q2a1c1b	Y5185
																								Q2a1c1b1	L301
																									Q2a1c1b1a	FGC4640
																						Q2a1c2	BZ3148
																							Q2a1c2a	BZ3149
																							Q2a1c2b	FGC27562
																			Q2b	L68.2
																				Q2b1	YP755
																					Q2b1a	YP754
																						Q2b1a1	Y15680
																							Q2b1a1b	BZ528
																				Q2b2	Y29468
																					Q2b2a	YP4500
																						Q2b2a1	Y28557
																							Q2b2a1a	Y28562
																								Q2b2a1a1	Y28555
																									Q2b2a1a1a	Y34111
																										Q2b2a1a1a1	BY185054
																										Q2b2a1a1a2	Y34109
																				Q2b3	Z36070
																	R	M207
																		R1	M173
																			R1a	L146
																				R1a1	M459
																					R1a1a	M512
																						R1a1a1	M417
																								R1a1a1a1	CTS7083
																							R1a1a1b	PF6162
																								R1a1a1b1	PF6217
																									R1a1a1b1a	S198
																											R1a1a1b1a1a	M458
																													R1a1a1b1a1a1a	L260
																														R1a1a1b1a1a1a1	YP256
																															R1a1a1b1a1a1a1a	YP254
																																R1a1a1b1a1a1a1a1	Y2905
																																R1a1a1b1a1a1a1a2	Y4135
																																R1a1a1b1a1a1a1a3	YP414
																													R1a1a1b1a1a1c	CTS11962.1
																														R1a1a1b1a1a1c1	L1029
																															R1a1a1b1a1a1c1a	YP263
																															R1a1a1b1a1a1c1b	YP416
																															R1a1a1b1a1a1c1c	FGC20517
																																		R1a1a1b1a1a1c1d1a1	YP445
																										R1a1a1b1a2	S466
																											R1a1a1b1a2a	S205
																															R1a1a1b1a2a1a1a	YP569
																												R1a1a1b1a2a2	YP270
																															R1a1a1b1a2a2a1a	YP350
																											R1a1a1b1a2b	CTS1211
																												R1a1a1b1a2b1	P278.2
																												R1a1a1b1a2b2	L784
																													R1a1a1b1a2b3a	CTS3402
																																		R1a1a1b1a2b3a1a1a1	L366
																																		R1a1a1b1a2b3a1a1a2	YP335
																																		R1a1a1b1a2b3a1a1b1	L365
																																						R1a1a1b1a2b3a1a1b1a1a1	F2686.1
																																	R1a1a1b1a2b3a3a1a	L1280
																																		R1a1a1b1a2b3a3a1a1	FGC11555
																																				R1a1a1b1a2b3a3a1b1a1	YP314
																																					R1a1a1b1a2b3a3a1b1a1a	YP331
																																	R1a1a1b1a2b3a3h	FGC10329
																											R1a1a1b1a2c	S24902
																											R1a1a1b1a3a	S221
																												R1a1a1b1a3a1	L448
																													R1a1a1b1a3a1a	CTS4179
																														R1a1a1b1a3a1a1	FGC11904
																															R1a1a1b1a3a1a1a	L176.1
																																R1a1a1b1a3a1a1a1	YP280
																																	R1a1a1b1a3a1a1a1a	FGC11917
																														R1a1a1b1a3a1a2	YP386
																															R1a1a1b1a3a1a3a	CTS3390
																												R1a1a1b1a3a2	S223
																													R1a1a1b1a3a2a	CTS8401
																														R1a1a1b1a3a2a1	CTS2243
																													R1a1a1b1a3a2b	CTS8277
																													R1a1a1b1a3a3a	CTS4027.1
																								R1a1a1b2	F992
																									R1a1a1b2a	F3105
																											R1a1a1b2a1a	L657.1
																													R1a1a1b2a1a1a	M787
																													R1a1a1b2a1a2c	AM00482
																										R1a1a1b2a2	Z2124
																												R1a1a1b2a2a1	Z2123
																											R1a1a1b2a2b	S4576
																												R1a1a1b2a2b1	F1345
																													R1a1a1b2a2b1a	CTS6
																													R1a1a1b2a2b1b	F2935
																												R1a1a1b2a2b2	Y57
																			R1b	M343
																				R1b1	L754
																					R1b1a	L388
																						R1b1a1	P297
																							R1b1a1a	M73
																								R1b1a1a1	M478
																									R1b1a1a1a	L1432
																									R1b1a1a1b	Y20750
																										R1b1a1a1b1	Y20747
																											R1b1a1a1b1a	Y22195
																								R1b1a1a2	BY15590
																							R1b1a1b	M269
																								R1b1a1b1	L23
																									R1b1a1b1a	L51
																										R1b1a1b1a1	P310
																											R1b1a1b1a1a	L151
																												R1b1a1b1a1a1	M405
																													R1b1a1b1a1a1a	FGC3861
																														R1b1a1b1a1a1a1	AM01881
																															R1b1a1b1a1a1a1a	L217.1
																															R1b1a1b1a1a1a1b	FGC17465
																																R1b1a1b1a1a1a1b1	FGC17464
																																	R1b1a1b1a1a1a1b1a	FGC17467
																														R1b1a1b1a1a1a2	FGC14875
																															R1b1a1b1a1a1a2a	FGC14874
																															R1b1a1b1a1a1a2b	A561
																																R1b1a1b1a1a1a2b1	FGC22182
																																R1b1a1b1a1a1a2b2	A560
																																	R1b1a1b1a1a1a2b2a	A565
																														R1b1a1b1a1a1a3	A1243
																													R1b1a1b1a1a1b	Z19
																														R1b1a1b1a1a1b1	Z14
																															R1b1a1b1a1a1b1a	S375
																																R1b1a1b1a1a1b1a1	L257
																																	R1b1a1b1a1a1b1a1a	S5741
																															R1b1a1b1a1a1b1b	DF95
																													R1b1a1b1a1a1c	S263
																														R1b1a1b1a1a1c1	S264
																															R1b1a1b1a1a1c1a	S497
																																R1b1a1b1a1a1c1a1	DF98
																																R1b1a1b1a1a1c1a2	DF96
																																	R1b1a1b1a1a1c1a2a	L1
																																		R1b1a1b1a1a1c1a2a1	A7108
																																		R1b1a1b1a1a1c1a2a2	A226
																																			R1b1a1b1a1a1c1a2a2a	A317
																																				R1b1a1b1a1a1c1a2a2a1	A9872
																																				R1b1a1b1a1a1c1a2a2a2	FGC20726
																																	R1b1a1b1a1a1c1a2b	FGC13326
																																		R1b1a1b1a1a1c1a2b1	S25234
																																			R1b1a1b1a1a1c1a2b1a	FGC13324
																																				R1b1a1b1a1a1c1a2b1a1	S11477
																																					R1b1a1b1a1a1c1a2b1a1a	P89.2
																																		R1b1a1b1a1a1c1a2b2	S22047
																																			R1b1a1b1a1a1c1a2b2a	FGC13602
																																				R1b1a1b1a1a1c1a2b2a1	FGC13595
																																					R1b1a1b1a1a1c1a2b2a1a	FGC13609
																															R1b1a1b1a1a1c1b	S5520
																																R1b1a1b1a1a1c1b1	S5556
																																R1b1a1b1a1a1c1b2	FGC11660
																														R1b1a1b1a1a1c2	S499
																															R1b1a1b1a1a1c2a	S1688
																																R1b1a1b1a1a1c2a1	M467
																																	R1b1a1b1a1a1c2a1a	Z37884
																																		R1b1a1b1a1a1c2a1a1	Z37885
																																	R1b1a1b1a1a1c2a1b	S16994
																																		R1b1a1b1a1a1c2a1b1	S16906
																																		R1b1a1b1a1a1c2a1b2	JFS0011
																																			R1b1a1b1a1a1c2a1b2a	JFS0010
																																	R1b1a1b1a1a1c2a1c	S15627
																																		R1b1a1b1a1a1c2a1c1	BY1053
																																			R1b1a1b1a1a1c2a1c1a	CTS4089.2
																																		R1b1a1b1a1a1c2a1c2	DF89
																																			R1b1a1b1a1a1c2a1c2a	FGC29371
																																				R1b1a1b1a1a1c2a1c2a1	JFS2001
																																			R1b1a1b1a1a1c2a1c2b	FGC12305
																																				R1b1a1b1a1a1c2a1c2b1	JFS0006
																																				R1b1a1b1a1a1c2a1c2b2	FGC12307
																																					R1b1a1b1a1a1c2a1c2b2a	FGC12312
																																					R1b1a1b1a1a1c2a1c2b2b	A6402
																																						R1b1a1b1a1a1c2a1c2b2b1	JFS0003
																																						R1b1a1b1a1a1c2a1c2b2b2	JFS0002
																																	R1b1a1b1a1a1c2a1d	DF93
																																		R1b1a1b1a1a1c2a1d1	JFS0009
																																			R1b1a1b1a1a1c2a1d1a	JFS0008
																																		R1b1a1b1a1a1c2a1d2	DF94
																																			R1b1a1b1a1a1c2a1d2a	S9787
																																				R1b1a1b1a1a1c2a1d2a1	JFS0007
																																				R1b1a1b1a1a1c2a1d2a2	S4060
																																				R1b1a1b1a1a1c2a1d2a3	S4078
																																					R1b1a1b1a1a1c2a1d2a3a	BY1351
																																					R1b1a1b1a1a1c2a1d2a3b	S4056
																																						R1b1a1b1a1a1c2a1d2a3b1	S4076
																																							R1b1a1b1a1a1c2a1d2a3b1a	S4057
																																R1b1a1b1a1a1c2a2	FGC19753
																															R1b1a1b1a1a1c2b	L48
																																R1b1a1b1a1a1c2b1	L47
																																	R1b1a1b1a1a1c2b1a	L44
																																		R1b1a1b1a1a1c2b1a1	L163
																																			R1b1a1b1a1a1c2b1a1a	L46
																																				R1b1a1b1a1a1c2b1a1a1	L45
																																					R1b1a1b1a1a1c2b1a1a1a	L477
																																						R1b1a1b1a1a1c2b1a1a1a1	FGC30616
																																					R1b1a1b1a1a1c2b1a1a1b	L292
																																			R1b1a1b1a1a1c2b1a1b	FGC6202
																																				R1b1a1b1a1a1c2b1a1b1	FGC6212
																																				R1b1a1b1a1a1c2b1a1b2	A6707
																																	R1b1a1b1a1a1c2b1b	Z159
																																		R1b1a1b1a1a1c2b1b1	S8368
																																		R1b1a1b1a1a1c2b1b2	FGC15332
																																			R1b1a1b1a1a1c2b1b2a	A1142
																																				R1b1a1b1a1a1c2b1b2a1	FGC15333
																																		R1b1a1b1a1a1c2b1b3	S9257
																																			R1b1a1b1a1a1c2b1b3a	S1260
																																			R1b1a1b1a1a1c2b1b3b	CTS6353
																																		R1b1a1b1a1a1c2b1b4	S3251.2
																																			R1b1a1b1a1a1c2b1b4a	FGC17296
																																				R1b1a1b1a1a1c2b1b4a1	FGC17308
																																				R1b1a1b1a1a1c2b1b4a2	FGC17297
																																					R1b1a1b1a1a1c2b1b4a2a	Y6456
																																						R1b1a1b1a1a1c2b1b4a3a2	FGC51456
																																			R1b1a1b1a1a1c2b1b4b	S6915
																																				R1b1a1b1a1a1c2b1b4b1	PH2129
																																				R1b1a1b1a1a1c2b1b4b2	FGC33315
																																				R1b1a1b1a1a1c2b1b4b3	S6925
																																			R1b1a1b1a1a1c2b1b4c	S20039
																																			R1b1a1b1a1a1c2b1b4d	A688
																																				R1b1a1b1a1a1c2b1b4d1	A6715
																																					R1b1a1b1a1a1c2b1b4d1a	CTS13009.2
																																					R1b1a1b1a1a1c2b1b4d1b	CTS3553
																																			R1b1a1b1a1a1c2b1b4e	FGC8590
																																				R1b1a1b1a1a1c2b1b4e1	FGC8579
																																					R1b1a1b1a1a1c2b1b4e1a	FGC8578
																																						R1b1a1b1a1a1c2b1b4e1a1	FGC13260
																																						R1b1a1b1a1a1c2b1b4e1a2	FGC8573
																																							R1b1a1b1a1a1c2b1b4e1a2a	A689
																																							R1b1a1b1a1a1c2b1b4e1a2b	S21495
																																							R1b1a1b1a1a1c2b1b4e1a2c	FGC8587
																																								R1b1a1b1a1a1c2b1b4e1a2c1	S378
																																R1b1a1b1a1a1c2b2	S268
																																	R1b1a1b1a1a1c2b2a	S271
																																		R1b1a1b1a1a1c2b2a1	Z32
																																			R1b1a1b1a1a1c2b2a1a	L844.2
																																			R1b1a1b1a1a1c2b2a1b	Z2
																																				R1b1a1b1a1a1c2b2a1b1	S272
																																					R1b1a1b1a1a1c2b2a1b1a	S515
																																						R1b1a1b1a1a1c2b2a1b1a1	S274
																																							R1b1a1b1a1a1c2b2a1b1a1a	S385
																																								R1b1a1b1a1a1c2b2a1b1a1a1	A5616
																																								R1b1a1b1a1a1c2b2a1b1a1a2	Z8175
																																									R1b1a1b1a1a1c2b2a1b1a1a2a	CTS10742
																																										R1b1a1b1a1a1c2b2a1b1a1a2a1	A574
																																									R1b1a1b1a1a1c2b2a1b1a1a2b	FGC12057
																																										R1b1a1b1a1a1c2b2a1b1a1a2b1	FGC29368
																																											R1b1a1b1a1a1c2b2a1b1a1a2b1a	L148
																																										R1b1a1b1a1a1c2b2a1b1a1a2b2	S18890
																																											R1b1a1b1a1a1c2b2a1b1a1a2b2a	FGC12058
																																												R1b1a1b1a1a1c2b2a1b1a1a2b2a1	A687
																																							R1b1a1b1a1a1c2b2a1b1a1b	S1774
																																								R1b1a1b1a1a1c2b2a1b1a1b1	M10044
																																								R1b1a1b1a1a1c2b2a1b1a1b2	S16701
																																									R1b1a1b1a1a1c2b2a1b1a1b2a	Y19921
																																									R1b1a1b1a1a1c2b2a1b1a1b2b	Z25289.2
																																						R1b1a1b1a1a1c2b2a1b1a2	FGC17517
																																						R1b1a1b1a1a1c2b2a1b1a3	S12035
																																							R1b1a1b1a1a1c2b2a1b1a3a	S9565
																																								R1b1a1b1a1a1c2b2a1b1a3a1	S26379
																																						R1b1a1b1a1a1c2b2a1b1a4	S275
																																							R1b1a1b1a1a1c2b2a1b1a4a	S514
																																								R1b1a1b1a1a1c2b2a1b1a4a1	S276
																																									R1b1a1b1a1a1c2b2a1b1a4a1a	S18951
																																										R1b1a1b1a1a1c2b2a1b1a4a1a1	CTS6428
																																											R1b1a1b1a1a1c2b2a1b1a4a1a1a	S10415
																																									R1b1a1b1a1a1c2b2a1b1a4a1b	A96
																																								R1b1a1b1a1a1c2b2a1b1a4a2	FGC29397
																																								R1b1a1b1a1a1c2b2a1b1a4a3	A300
																																							R1b1a1b1a1a1c2b2a1b1a4b	S512
																																								R1b1a1b1a1a1c2b2a1b1a4b1	S387
																																									R1b1a1b1a1a1c2b2a1b1a4b1a	CTS5601
																																										R1b1a1b1a1a1c2b2a1b1a4b1a1	FGC4150
																																											R1b1a1b1a1a1c2b2a1b1a4b1a1a	CTS7080
																																												R1b1a1b1a1a1c2b2a1b1a4b1a1a1	FGC4161
																																												R1b1a1b1a1a1c2b2a1b1a4b1a1a2	S14406
																																													R1b1a1b1a1a1c2b2a1b1a4b1a1a2a	Y15631
																																														R1b1a1b1a1a1c2b2a1b1a4b1a1a2a1	Y15633
																																													R1b1a1b1a1a1c2b2a1b1a4b1a1a2b	BY118
																																									R1b1a1b1a1a1c2b2a1b1a4b1b	FGC36939
																																									R1b1a1b1a1a1c2b2a1b1a4b1c	FGC11784
																																										R1b1a1b1a1a1c2b2a1b1a4b1c1	S6881
																																											R1b1a1b1a1a1c2b2a1b1a4b1c1a	A8050
																																								R1b1a1b1a1a1c2b2a1b1a4b2	DF101
																																									R1b1a1b1a1a1c2b2a1b1a4b2a	DF102
																																										R1b1a1b1a1a1c2b2a1b1a4b2a1	A6906
																																										R1b1a1b1a1a1c2b2a1b1a4b2a2	S5245
																																											R1b1a1b1a1a1c2b2a1b1a4b2a2a	Z27210
																																											R1b1a1b1a1a1c2b2a1b1a4b2a2b	FGC12993
																																												R1b1a1b1a1a1c2b2a1b1a4b2a2b1	A321
																																											R1b1a1b1a1a1c2b2a1b1a4b2a2c	S5627
																																												R1b1a1b1a1a1c2b2a1b1a4b2a2c1	S5246
																																													R1b1a1b1a1a1c2b2a1b1a4b2a2c1a	FGC15254
																																														R1b1a1b1a1a1c2b2a1b1a4b2a2c1a1	FGC35613
																																						R1b1a1b1a1a1c2b2a1b1a5	FGC20422
																																						R1b1a1b1a1a1c2b2a1b1a6	FGC5264
																																							R1b1a1b1a1a1c2b2a1b1a6a	FGC5254
																																								R1b1a1b1a1a1c2b2a1b1a6a1	FGC5253
																																					R1b1a1b1a1a1c2b2a1b1b	FGC7559
																																						R1b1a1b1a1a1c2b2a1b1b1	CTS10893
																																							R1b1a1b1a1a1c2b2a1b1b1a	FGC924
																																								R1b1a1b1a1a1c2b2a1b1b1a1	FGC918
																																									R1b1a1b1a1a1c2b2a1b1b1a1a	Y10790
																																									R1b1a1b1a1a1c2b2a1b1b1a1b	FGC934
																																										R1b1a1b1a1a1c2b2a1b1b1a1b1	A6896
																																										R1b1a1b1a1a1c2b2a1b1b1a1b2	FGC921
																																											R1b1a1b1a1a1c2b2a1b1b1a1b2a	FGC925
																																											R1b1a1b1a1a1c2b2a1b1b1a1b2b	FGC17852
																																								R1b1a1b1a1a1c2b2a1b1b1a2	A294
																																									R1b1a1b1a1a1c2b2a1b1b1a2a	Y10791
																																									R1b1a1b1a1a1c2b2a1b1b1a2b	A295
																																									R1b1a1b1a1a1c2b2a1b1b1a2c	A6886
																																								R1b1a1b1a1a1c2b2a1b1b1a3	A8151
																																									R1b1a1b1a1a1c2b2a1b1b1a3a	A9455
																																										R1b1a1b1a1a1c2b2a1b1b1a3a1	A9454
																																								R1b1a1b1a1a1c2b2a1b1b1a4	A6389
																																									R1b1a1b1a1a1c2b2a1b1b1a4a	BY3323
																																									R1b1a1b1a1a1c2b2a1b1b1a4b	S27458
																																								R1b1a1b1a1a1c2b2a1b1b1a5	FGC17429
																																									R1b1a1b1a1a1c2b2a1b1b1a5a	S6898
																																									R1b1a1b1a1a1c2b2a1b1b1a5b	FGC14717
																																					R1b1a1b1a1a1c2b2a1b1c	FGC17344
																																				R1b1a1b1a1a1c2b2a1b2	S22165
																																				R1b1a1b1a1a1c2b2a1b3	S16218
																																				R1b1a1b1a1a1c2b2a1b4	S20054
																																				R1b1a1b1a1a1c2b2a1b5	S15510
																																					R1b1a1b1a1a1c2b2a1b5a	Y7369
																																						R1b1a1b1a1a1c2b2a1b5a1	Y7370
																																						R1b1a1b1a1a1c2b2a1b5a2	S25738
																																				R1b1a1b1a1a1c2b2a1b6	S9342
																																					R1b1a1b1a1a1c2b2a1b6a	PF5143
																																						R1b1a1b1a1a1c2b2a1b6a1	PH765
																																					R1b1a1b1a1a1c2b2a1b6b	CTS9539
																																						R1b1a1b1a1a1c2b2a1b6b1	S10957
																																	R1b1a1b1a1a1c2b2b	S1743
																																		R1b1a1b1a1a1c2b2b1	S505
																																			R1b1a1b1a1a1c2b2b1a	Z326
																																				R1b1a1b1a1a1c2b2b1a1	FGC10367
																																					R1b1a1b1a1a1c2b2b1a1a	S1741
																																						R1b1a1b1a1a1c2b2b1a1a1	CTS2509
																																							R1b1a1b1a1a1c2b2b1a1a1a	BY556
																																							R1b1a1b1a1a1c2b2b1a1a1b	FGC363
																																								R1b1a1b1a1a1c2b2b1a1a1b1	BY557
																																								R1b1a1b1a1a1c2b2b1a1a1b2	S381
																																									R1b1a1b1a1a1c2b2b1a1a1b2a	L1251
																																									R1b1a1b1a1a1c2b2b1a1a1b2b	S508
																																							R1b1a1b1a1a1c2b2b1a1a1c	FGC13492
																																							R1b1a1b1a1a1c2b2b1a1a1d	FGC15317
																																								R1b1a1b1a1a1c2b2b1a1a1d1	FGC15315
																																								R1b1a1b1a1a1c2b2b1a1a1d2	A2269
																																							R1b1a1b1a1a1c2b2b1a1a1e	FGC564
																																								R1b1a1b1a1a1c2b2b1a1a1e1	S10817
																																								R1b1a1b1a1a1c2b2b1a1a1e2	CTS9
																																								R1b1a1b1a1a1c2b2b1a1a1e3	Z5054
																																									R1b1a1b1a1a1c2b2b1a1a1e3a	FGC539
																																										R1b1a1b1a1a1c2b2b1a1a1e3a1	L188
																																				R1b1a1b1a1a1c2b2b1a2	S21728
																																	R1b1a1b1a1a1c2b2c	S22294
																																	R1b1a1b1a1a1c2b2d	S4111
																																R1b1a1b1a1a1c2b3	S23189
																																	R1b1a1b1a1a1c2b3a	L200
																																		R1b1a1b1a1a1c2b3a1	S9355
																																			R1b1a1b1a1a1c2b3a1a	S21607
																																				R1b1a1b1a1a1c2b3a1a1	A5317
																																					R1b1a1b1a1a1c2b3a1a2a	A765
																																	R1b1a1b1a1a1c2b3b	S10271
																																	R1b1a1b1a1a1c2b3c	A6706
																																		R1b1a1b1a1a1c2b3c1	BY3320
																																		R1b1a1b1a1a1c2b3c2	S19342
																																			R1b1a1b1a1a1c2b3c2a	S18372
																															R1b1a1b1a1a1c2c	FGC20667
																																R1b1a1b1a1a1c2c1	FGC20676
																															R1b1a1b1a1a1c2d	FGC8512
																																R1b1a1b1a1a1c2d1	Z155
																																	R1b1a1b1a1a1c2d1a	FGC8507
																															R1b1a1b1a1a1c2e	FGC13938
																																R1b1a1b1a1a1c2e1	A7220
																																R1b1a1b1a1a1c2e2	S9891
																																	R1b1a1b1a1a1c2e2a	S10250
																														R1b1a1b1a1a1c3	M323.2
																													R1b1a1b1a1a1d	FGC396
																														R1b1a1b1a1a1d1	L199.1
																														R1b1a1b1a1a1d2	Y11145
																													R1b1a1b1a1a1e	S12025
																														R1b1a1b1a1a1e1	S16361
																															R1b1a1b1a1a1e1a	S25186
																															R1b1a1b1a1a1e1b	FGC15048
																														R1b1a1b1a1a1e2	S15526
																															R1b1a1b1a1a1e2a	S25007
																													R1b1a1b1a1a1f	A2150
																														R1b1a1b1a1a1f1	Y9130
																													R1b1a1b1a1a1g	S18632
																														R1b1a1b1a1a1g1	S21017
																														R1b1a1b1a1a1g2	FGC35808
																													R1b1a1b1a1a1h	S11493
																												R1b1a1b1a1a2	P312
																													R1b1a1b1a1a2a	DF27
																														R1b1a1b1a1a2a1	Z195
																															R1b1a1b1a1a2a1a	Z272
																																R1b1a1b1a1a2a1a1	S230
																																	R1b1a1b1a1a2a1a1a	S1217
																																		R1b1a1b1a1a2a1a1a1	S181
																																			R1b1a1b1a1a2a1a1a1a	S348
																																				R1b1a1b1a1a2a1a1a1a1	M153
																																		R1b1a1b1a1a2a1a1a2	CTS4065
																																		R1b1a1b1a1a2a1a1a3	CTS3702
																																			R1b1a1b1a1a2a1a1a3a	BY42727
																																				R1b1a1b1a1a2a1a1a3a1	BY107537
																																					R1b1a1b1a1a2a1a1a3a1a	FT93521
																																R1b1a1b1a1a2a1a2	DF17
																															R1b1a1b1a1a2a1b	S228
																																R1b1a1b1a1a2a1b1	S460
																																	R1b1a1b1a1a2a1b1a	Z262
																																		R1b1a1b1a1a2a1b1a1	M167
																																R1b1a1b1a1a2a1b2	L165
																																R1b1a1b1a1a2a1b3	CTS4188
																														R1b1a1b1a1a2a2	Z2552
																															R1b1a1b1a1a2a2a	L617
																														R1b1a1b1a1a2a3	L881
																														R1b1a1b1a1a2a4	A431
																														R1b1a1b1a1a2a5	F1343
																														R1b1a1b1a1a2a6	Z2571
																													R1b1a1b1a1a2b	PF6570
																														R1b1a1b1a1a2b1	L2
																															R1b1a1b1a1a2b1a	S255
																																R1b1a1b1a1a2b1a1	L20
																																R1b1a1b1a1a2b1a2	S368
																																	R1b1a1b1a1a2b1a2a	S487
																																		R1b1a1b1a1a2b1a2a1	Z275
																															R1b1a1b1a1a2b1b	L196
																															R1b1a1b1a1a2b1c	Z49
																																R1b1a1b1a1a2b1c1	S211
																																	R1b1a1b1a1a2b1c1a	S369
																																		R1b1a1b1a1a2b1c1a1	L562
																																			R1b1a1b1a1a2b1c1a1a	S1468
																																				R1b1a1b1a1a2b1c1a1a1	Z147
																																					R1b1a1b1a1a2b1c1a1a1a	Z52
																																						R1b1a1b1a1a2b1c1a1a1a1	CTS11232
																																						R1b1a1b1a1a2b1c1a1a1a2	CTS278
																																							R1b1a1b1a1a2b1c1a1a1a2a	CTS4521
																																								R1b1a1b1a1a2b1c1a1a1a2a1	F3799.2
																																				R1b1a1b1a1a2b1c1a1a2	S1491
																																					R1b1a1b1a1a2b1c1a1a2a	S1480
																																						R1b1a1b1a1a2b1c1a1a2a1	Z29987
																																							R1b1a1b1a1a2b1c1a1a2a1a	Z29988
																																								R1b1a1b1a1a2b1c1a1a2a1a1	Z29989
																																			R1b1a1b1a1a2b1c1a1b	F1947.2
																																			R1b1a1b1a1a2b1c1a1c	PF27
																																		R1b1a1b1a1a2b1c1a2	F1036.2
																																	R1b1a1b1a1a2b1c1b	S257
																																		R1b1a1b1a1a2b1c1b1	CTS9490
																																			R1b1a1b1a1a2b1c1b1a	CTS7970
																																				R1b1a1b1a1a2b1c1b1a1	S7409
																																					R1b1a1b1a1a2b1c1b1a1a	FGC5007
																																				R1b1a1b1a1a2b1c1b1a2	CTS8125
																																					R1b1a1b1a1a2b1c1b1a2a	CTS11381
																																		R1b1a1b1a1a2b1c1b2	L654.1
																																			R1b1a1b1a1a2b1c1b2a	S42
																																				R1b1a1b1a1a2b1c1b2a1	CTS7197
																																		R1b1a1b1a1a2b1c1b3	FGC12378
																																			R1b1a1b1a1a2b1c1b3a	FGC12383
																																				R1b1a1b1a1a2b1c1b3a1	FGC12401
																																		R1b1a1b1a1a2b1c1b4	L553
																																		R1b1a1b1a1a2b1c1b5	L552
																																		R1b1a1b1a1a2b1c1b6	CTS1789
																																	R1b1a1b1a1a2b1c1c	FGC22940
																																		R1b1a1b1a1a2b1c1c1	FGC22952
																																R1b1a1b1a1a2b1c2	S8183
																																	R1b1a1b1a1a2b1c2a	S22778
																																		R1b1a1b1a1a2b1c2a1	S20782
																																	R1b1a1b1a1a2b1c2b	FGC31474
																																		R1b1a1b1a1a2b1c2b1	FGC20796
																																			R1b1a1b1a1a2b1c2b1a	A1168
																																				R1b1a1b1a1a2b1c2b1a1	FGC20777
																																					R1b1a1b1a1a2b1c2b1a1a	FGC20775
																																				R1b1a1b1a1a2b1c2b1a2	S8172
																																		R1b1a1b1a1a2b1c2b2	Y11179
																																			R1b1a1b1a1a2b1c2b2a	Y12696
																																				R1b1a1b1a1a2b1c2b2a1	Y13610
																															R1b1a1b1a1a2b1d	FGC22501
																																R1b1a1b1a1a2b1d1	FGC22538
																																		R1b1a1b1a1a2b1d1b1	FGC42117
																																			R1b1a1b1a1a2b1d1a3	FGC42109
																														R1b1a1b1a1a2b2	S206
																														R1b1a1b1a1a2b3	PF6601
																															R1b1a1b1a1a2b3a	L4
																															R1b1a1b1a1a2b3b	S47
																															R1b1a1b1a1a2b3c	PF6578
																													R1b1a1b1a1a2c	S461
																														R1b1a1b1a1a2c1	L21
																															R1b1a1b1a1a2c1a	CTS241
																																R1b1a1b1a1a2c1a1	Z39589
																																	R1b1a1b1a1a2c1a1a	DF49
																																		R1b1a1b1a1a2c1a1a1	S6154
																																			R1b1a1b1a1a2c1a1a1a	S476
																																				R1b1a1b1a1a2c1a1a1a1	DF23
																																					R1b1a1b1a1a2c1a1a1a1a	Z2961
																																						R1b1a1b1a1a2c1a1a1a1a1	M222
																																							R1b1a1b1a1a2c1a1a1a1a1a	DF106
																																								R1b1a1b1a1a2c1a1a1a1a1a1	DF104
																																									R1b1a1b1a1a2c1a1a1a1a1a1a	DF109
																																										R1b1a1b1a1a2c1a1a1a1a1a1a1	DF85
																																											R1b1a1b1a1a2c1a1a1a1a1a1a1a	S673
																																												R1b1a1b1a1a2c1a1a1a1a1a1a1a1	S668
																																													R1b1a1b1a1a2c1a1a1a1a1a1a1a1a	DF97
																																														R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1	FGC19851
																																															R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a	FGC19846
																																																R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a1	FGC19844
																																																	R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a1a	FGC19845
																																																		R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a1a1	A1147
																																																		R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a1a2	FGC19842
																																																		R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a1a3	FGC30115
																																																R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1a2	FGC42106
																																															R1b1a1b1a1a2c1a1a1a1a1a1a1a1a1b	Z29319
																																														R1b1a1b1a1a2c1a1a1a1a1a1a1a1a2	A1330
																																															R1b1a1b1a1a2c1a1a1a1a1a1a1a1a2a	A1332
																																													R1b1a1b1a1a2c1a1a1a1a1a1a1a2a	FGC32825
																																											R1b1a1b1a1a2c1a1a1a1a1a1a1b	FGC12183
																																											R1b1a1b1a1a2c1a1a1a1a1a1a1d	FGC35551
																																												R1b1a1b1a1a2c1a1a1a1a1a1a1d1	FGC35553
																																										R1b1a1b1a1a2c1a1a1a1a1a1a2	FGC4108
																																											R1b1a1b1a1a2c1a1a1a1a1a1a2a	S7814
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2a1	F1265.2
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2a2	A694
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2a3	BY205
																																													R1b1a1b1a1a2c1a1a1a1a1a1a2a3a	BY4142
																																											R1b1a1b1a1a2c1a1a1a1a1a1a2b	FGC7928
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2b1	FGC23592
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2b2	FGC4113
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2b3	FGC19860
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2b4	S590
																																													R1b1a1b1a1a2c1a1a1a1a1a1a2b4a	S595
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2b5	A9885
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2b6	A8591
																																											R1b1a1b1a1a2c1a1a1a1a1a1a2c	FGC30690
																																												R1b1a1b1a1a2c1a1a1a1a1a1a2c1	A10535
																																													R1b1a1b1a1a2c1a1a1a1a1a1a2c1a	FGC30697
																																										R1b1a1b1a1a2c1a1a1a1a1a1a3	A259
																																											R1b1a1b1a1a2c1a1a1a1a1a1a3a	A260
																																										R1b1a1b1a1a2c1a1a1a1a1a1a4	A223
																																											R1b1a1b1a1a2c1a1a1a1a1a1a4a	A224
																																											R1b1a1b1a1a2c1a1a1a1a1a1a4b	A822
																																												R1b1a1b1a1a2c1a1a1a1a1a1a4b1	A984
																																													R1b1a1b1a1a2c1a1a1a1a1a1a4b1a	A11242
																																										R1b1a1b1a1a2c1a1a1a1a1a1a5	FGC4133
																																							R1b1a1b1a1a2c1a1a1a1a1b	FGC4077
																																								R1b1a1b1a1a2c1a1a1a1a1b1	A725
																																								R1b1a1b1a1a2c1a1a1a1a1b2	FGC12948
																																							R1b1a1b1a1a2c1a1a1a1a1c	FGC451
																																								R1b1a1b1a1a2c1a1a1a1a1c1	FGC453
																																			R1b1a1b1a1a2c1a1a1b	FGC9631
																																				R1b1a1b1a1a2c1a1a1b1	ZP41
																																				R1b1a1b1a1a2c1a1a1b2	FGC9645
																																		R1b1a1b1a1a2c1a1a2	ZP20
																																			R1b1a1b1a1a2c1a1a2a	ZP21
																																				R1b1a1b1a1a2c1a1a2a1	ZP23
																																	R1b1a1b1a1a2c1a1b	CTS1751
																																	R1b1a1b1a1a2c1a1c	L371
																																	R1b1a1b1a1a2c1a1d	CTS2501
																																		R1b1a1b1a1a2c1a1d1	S775
																																			R1b1a1b1a1a2c1a1d1a	L744
																																				R1b1a1b1a1a2c1a1d1a1	S781
																																				R1b1a1b1a1a2c1a1d1a2	Y138948
																																				R1b1a1b1a1a2c1a1d1a3	Z38845
																																			R1b1a1b1a1a2c1a1d1b	A600
																																		R1b1a1b1a1a2c1a1d2	BY114
																																			R1b1a1b1a1a2c1a1d2a	L563
																																		R1b1a1b1a1a2c1a1d3	FGC5572
																																			R1b1a1b1a1a2c1a1d3a	MC21
																																				R1b1a1b1a1a2c1a1d3a1	FGC36636
																																				R1b1a1b1a1a2c1a1d3a2	BY39380
																																			R1b1a1b1a1a2c1a1d3b	FGC5586
																																				R1b1a1b1a1a2c1a1d3b1	A40
																																					R1b1a1b1a1a2c1a1d3b1a	FGC5584
																																						R1b1a1b1a1a2c1a1d3b1a1	FGC5581
																																					R1b1a1b1a1a2c1a1d3b1b	MC04
																																						R1b1a1b1a1a2c1a1d3b1b1	MC13
																																						R1b1a1b1a1a2c1a1d3b1b2	MC15
																																					R1b1a1b1a1a2c1a1d3b1c	BY281
																																			R1b1a1b1a1a2c1a1d3c	BY166
																																			R1b1a1b1a1a2c1a1d3d	BY23462
																																		R1b1a1b1a1a2c1a1d4	A874
																																			R1b1a1b1a1a2c1a1d4a	A875
																																				R1b1a1b1a1a2c1a1d4a1	A1786
																																					R1b1a1b1a1a2c1a1d4a1a	F1680
																																		R1b1a1b1a1a2c1a1d5	A98
																																			R1b1a1b1a1a2c1a1d5a	A97
																																		R1b1a1b1a1a2c1a1d6	FGC13023
																																			R1b1a1b1a1a2c1a1d6a	FGC13017
																																				R1b1a1b1a1a2c1a1d6a1	FGC13014
																																		R1b1a1b1a1a2c1a1d7	Y5628
																																	R1b1a1b1a1a2c1a1e	S470
																																		R1b1a1b1a1a2c1a1e1	L555
																																	R1b1a1b1a1a2c1a1f	L1335
																																		R1b1a1b1a1a2c1a1f1	CTS11722
																																			R1b1a1b1a1a2c1a1f1a	S744
																																				R1b1a1b1a1a2c1a1f1a1	S764
																																					R1b1a1b1a1a2c1a1f1a1a	L743
																																					R1b1a1b1a1a2c1a1f1a1b	FGC17596
																																						R1b1a1b1a1a2c1a1f1a1b1	FGC17595
																																					R1b1a1b1a1a2c1a1f1a1c	S756
																																				R1b1a1b1a1a2c1a1f1a2	S691
																																				R1b1a1b1a1a2c1a1f1a3	PF5236
																																				R1b1a1b1a1a2c1a1f1a4	CTS4931
																																			R1b1a1b1a1a2c1a1f1b	Z16329
																																			R1b1a1b1a1a2c1a1f1c	FGC10125
																																				R1b1a1b1a1a2c1a1f1c1	FGC10117
																																	R1b1a1b1a1a2c1a1g	S16264
																																		R1b1a1b1a1a2c1a1g1	L679
																																	R1b1a1b1a1a2c1a1h	Y5305
																																		R1b1a1b1a1a2c1a1h1	MC02
																																			R1b1a1b1a1a2c1a1h1a	A7309
																																	R1b1a1b1a1a2c1a1i	S1051
																																		R1b1a1b1a1a2c1a1i1	FGC17938
																																		R1b1a1b1a1a2c1a1i2	S1050
																																		R1b1a1b1a1a2c1a1i3	FGC17906
																																			R1b1a1b1a1a2c1a1i3a	FGC17907
																																			R1b1a1b1a1a2c1a1i3b	FGC29039
																																	R1b1a1b1a1a2c1a1j	FGC13780
																																	R1b1a1b1a1a2c1a1k	S1026
																																	R1b1a1b1a1a2c1a1l	FGC35995
																																	R1b1a1b1a1a2c1a1m	A4556
																																	R1b1a1b1a1a2c1a1n	BY575
																																R1b1a1b1a1a2c1a2	CTS5396
																																	R1b1a1b1a1a2c1a2a	S5668
																																		R1b1a1b1a1a2c1a2a1	Z16340
																																			R1b1a1b1a1a2c1a2a1a	FGC9807
																																				R1b1a1b1a1a2c1a2a1a1	FGC9793
																																					R1b1a1b1a1a2c1a2a1a1a	FGC9804
																																						R1b1a1b1a1a2c1a2a1a1a1	FGC9800
																																							R1b1a1b1a1a2c1a2a1a1a1a	P66_1
																																							R1b1a1b1a1a2c1a2a1a1a1b	FGC9808
																																								R1b1a1b1a1a2c1a2a1a1a1b1	BY2598
																																						R1b1a1b1a1a2c1a2a1a1a2	FGC32571
																																					R1b1a1b1a1a2c1a2a1a1b	Z20655
																																					R1b1a1b1a1a2c1a2a1a1c	ZS347
																																		R1b1a1b1a1a2c1a2a2	A7
																																			R1b1a1b1a1a2c1a2a2a	S5979
																																				R1b1a1b1a1a2c1a2a2a1	S5982
																																					R1b1a1b1a1a2c1a2a2a1a	A3
																																						R1b1a1b1a1a2c1a2a2a1a1	A8
																																							R1b1a1b1a1a2c1a2a2a1a1a	Y11636
																																						R1b1a1b1a1a2c1a2a2a1a2	Z18065
																																					R1b1a1b1a1a2c1a2a2a1b	Z17816
																																						R1b1a1b1a1a2c1a2a2a1b1	Z17815
																																							R1b1a1b1a1a2c1a2a2a1b1a	A1067
																																							R1b1a1b1a1a2c1a2a2a1b1b	A1075
																																					R1b1a1b1a1a2c1a2a2a1c	Z17817
																																						R1b1a1b1a1a2c1a2a2a1c1	BY615
																																							R1b1a1b1a1a2c1a2a2a1c1a	BY451
																																					R1b1a1b1a1a2c1a2a2a1d	BY2570
																																						R1b1a1b1a1a2c1a2a2a1d1	BY651
																																							R1b1a1b1a1a2c1a2a2a1d1a	BY2634
																																					R1b1a1b1a1a2c1a2a2a1e	ZS4585
																																			R1b1a1b1a1a2c1a2a2b	Z17819
																																				R1b1a1b1a1a2c1a2a2b1	FGC13508
																																					R1b1a1b1a1a2c1a2a2b1a	Z18489
																																						R1b1a1b1a1a2c1a2a2b1a1	Z18490
																																							R1b1a1b1a1a2c1a2a2b1a1a	BY2630
																																			R1b1a1b1a1a2c1a2a2c	Z17833
																																			R1b1a1b1a1a2c1a2a2d	S7834
																																				R1b1a1b1a1a2c1a2a2d1	Z21262
																																					R1b1a1b1a1a2c1a2a2d1a	Z21269
																																						R1b1a1b1a1a2c1a2a2d1a1	L577
																																		R1b1a1b1a1a2c1a2a3	Z16343
																																			R1b1a1b1a1a2c1a2a3a	Z16855
																																			R1b1a1b1a1a2c1a2a3b	Z17911
																																	R1b1a1b1a1a2c1a2b	S6365
																																		R1b1a1b1a1a2c1a2b1	CTS6621
																																			R1b1a1b1a1a2c1a2b1a	Z16403
																																				R1b1a1b1a1a2c1a2b1a1	Z16407
																																					R1b1a1b1a1a2c1a2b1a1a	Z17909
																																						R1b1a1b1a1a2c1a2b1a1a1	BY2745
																																						R1b1a1b1a1a2c1a2b1a1a2	ZW04
																																							R1b1a1b1a1a2c1a2b1a1a2a	ZW02
																																						R1b1a1b1a1a2c1a2b1a1a3	BY55195
																																					R1b1a1b1a1a2c1a2b1a1b	BY553
																																		R1b1a1b1a1a2c1a2b2	BY16
																																			R1b1a1b1a1a2c1a2b2a	CTS3087
																																				R1b1a1b1a1a2c1a2b2a1	Z17626
																																					R1b1a1b1a1a2c1a2b2a1a	Z17624
																																			R1b1a1b1a1a2c1a2b2b	Z16372
																																				R1b1a1b1a1a2c1a2b2b1	BY182
																																					R1b1a1b1a1a2c1a2b2b1a	M2775
																																				R1b1a1b1a1a2c1a2b2b2	BY404
																																		R1b1a1b1a1a2c1a2b3	BY17
																																			R1b1a1b1a1a2c1a2b3a	FGC11794
																																				R1b1a1b1a1a2c1a2b3a1	FGC11801
																																					R1b1a1b1a1a2c1a2b3a1a	FGC11788
																																	R1b1a1b1a1a2c1a2c	Z23534
																																	R1b1a1b1a1a2c1a2d	FGC13437
																																R1b1a1b1a1a2c1a3	FGC11134
																																	R1b1a1b1a1a2c1a3a	A146
																																		R1b1a1b1a1a2c1a3a1	L96
																																		R1b1a1b1a1a2c1a3a2	CTS4466
																																			R1b1a1b1a1a2c1a3a2a	S1115
																																				R1b1a1b1a1a2c1a3a2a1	A541
																																					R1b1a1b1a1a2c1a3a2a1a	S1121
																																						R1b1a1b1a1a2c1a3a2a1a1	Z16521
																																							R1b1a1b1a1a2c1a3a2a1a1a	A1133
																																						R1b1a1b1a1a2c1a3a2a1a2	FGC11147
																																							R1b1a1b1a1a2c1a3a2a1a2a	FGC29280
																																								R1b1a1b1a1a2c1a3a2a1a2a1	S1126
																																									R1b1a1b1a1a2c1a3a2a1a2a1a	A542
																																									R1b1a1b1a1a2c1a3a2a1a2a1b	A726
																																								R1b1a1b1a1a2c1a3a2a1a2a2	FGC29290
																																							R1b1a1b1a1a2c1a3a2a1a2b	A150
																																								R1b1a1b1a1a2c1a3a2a1a2b1	FGC11145
																																								R1b1a1b1a1a2c1a3a2a1a2b2	A745
																																							R1b1a1b1a1a2c1a3a2a1a2c	A159
																																								R1b1a1b1a1a2c1a3a2a1a2c1	BY149
																																									R1b1a1b1a1a2c1a3a2a1a2c1a	A664
																																									R1b1a1b1a1a2c1a3a2a1a2c1b	A923
																																									R1b1a1b1a1a2c1a3a2a1a2c1c	A6506
																																								R1b1a1b1a1a2c1a3a2a1a2c2	A9185.2
																																							R1b1a1b1a1a2c1a3a2a1a2d	BY2880
																																								R1b1a1b1a1a2c1a3a2a1a2d1	A804
																																									R1b1a1b1a1a2c1a3a2a1a2d1a	A802
																																									R1b1a1b1a1a2c1a3a2a1a2d1b	FGC29067
																																										R1b1a1b1a1a2c1a3a2a1a2d1b1	FGC29071
																																								R1b1a1b1a1a2c1a3a2a1a2d2	A2221
																																									R1b1a1b1a1a2c1a3a2a1a2d2a	A6464
																																							R1b1a1b1a1a2c1a3a2a1a2e	A806
																																					R1b1a1b1a1a2c1a3a2a1b	A195
																																						R1b1a1b1a1a2c1a3a2a1b1	A761
																																							R1b1a1b1a1a2c1a3a2a1b1a	A88
																																								R1b1a1b1a1a2c1a3a2a1b1a1	Z16259
																																									R1b1a1b1a1a2c1a3a2a1b1a1a	ZS4589
																																									R1b1a1b1a1a2c1a3a2a1b1a1b	A1333
																																									R1b1a1b1a1a2c1a3a2a1b1a1c	A2220
																																										R1b1a1b1a1a2c1a3a2a1b1a1c1	A2219
																																											R1b1a1b1a1a2c1a3a2a1b1a1c1a	A2294
																																								R1b1a1b1a1a2c1a3a2a1b1a2	A89
																																									R1b1a1b1a1a2c1a3a2a1b1a2a	A155.2
																																										R1b1a1b1a1a2c1a3a2a1b1a2a1	A156
																																								R1b1a1b1a1a2c1a3a2a1b1a3	A7660
																																							R1b1a1b1a1a2c1a3a2a1b1b	Z16254
																																								R1b1a1b1a1a2c1a3a2a1b1b1	A154
																																									R1b1a1b1a1a2c1a3a2a1b1b1a	A153
																																										R1b1a1b1a1a2c1a3a2a1b1b1a1	A474
																																					R1b1a1b1a1a2c1a3a2a1c	A956
																																						R1b1a1b1a1a2c1a3a2a1c1	A2288
																																							R1b1a1b1a1a2c1a3a2a1c1a	A5312
																																								R1b1a1b1a1a2c1a3a2a1c1a1	A7753
																																						R1b1a1b1a1a2c1a3a2a1c2	A7754
																																					R1b1a1b1a1a2c1a3a2a1d	A151
																																						R1b1a1b1a1a2c1a3a2a1d1	A661
																																							R1b1a1b1a1a2c1a3a2a1d1a	BY140
																																						R1b1a1b1a1a2c1a3a2a1d2	A714
																																						R1b1a1b1a1a2c1a3a2a1d3	B42
																																				R1b1a1b1a1a2c1a3a2a2	A212
																																					R1b1a1b1a1a2c1a3a2a2a	A206
																																					R1b1a1b1a1a2c1a3a2a2b	A7699
																																						R1b1a1b1a1a2c1a3a2a2b1	A6525
																																							R1b1a1b1a1a2c1a3a2a2b1a	A7756
																																						R1b1a1b1a1a2c1a3a2a2b2	A6526
																																					R1b1a1b1a1a2c1a3a2a2c	A7648
																																				R1b1a1b1a1a2c1a3a2a3	A663
																																					R1b1a1b1a1a2c1a3a2a3a	A2289
																																			R1b1a1b1a1a2c1a3a2b	A7751
																																	R1b1a1b1a1a2c1a3b	A286
																																		R1b1a1b1a1a2c1a3b1	FGC18030
																																	R1b1a1b1a1a2c1a3c	FGC11293
																																R1b1a1b1a1a2c1a4	ZZ10
																																	R1b1a1b1a1a2c1a4a	S219
																																		R1b1a1b1a1a2c1a4a1	L159.2
																																			R1b1a1b1a1a2c1a4a1a	Z16433
																																				R1b1a1b1a1a2c1a4a1a1	A1248
																																				R1b1a1b1a1a2c1a4a1a2	A430
																																					R1b1a1b1a1a2c1a4a1a2a	A6308
																																						R1b1a1b1a1a2c1a4a1a2a1	A476
																																					R1b1a1b1a1a2c1a4a1a2b	A1367
																																						R1b1a1b1a1a2c1a4a1a2b1	FGC32914
																																							R1b1a1b1a1a2c1a4a1a2b1a	A816
																																						R1b1a1b1a1a2c1a4a1a2b2	A5311
																																					R1b1a1b1a1a2c1a4a1a2c	A1154
																																						R1b1a1b1a1a2c1a4a1a2c1	A8430
																																	R1b1a1b1a1a2c1a4b	S218
																																		R1b1a1b1a1a2c1a4b1	L554
																																		R1b1a1b1a1a2c1a4b2	S868
																																			R1b1a1b1a1a2c1a4b2a	L226
																																				R1b1a1b1a1a2c1a4b2a1	FGC5660
																																					R1b1a1b1a1a2c1a4b2a1a	FGC5628
																																						R1b1a1b1a1a2c1a4b2a1a1	DC1
																																							R1b1a1b1a1a2c1a4b2a1a1a	Y6913
																																						R1b1a1b1a1a2c1a4b2a1a2	DC24
																																					R1b1a1b1a1a2c1a4b2a1b	FGC12290
																																					R1b1a1b1a1a2c1a4b2a1c	DC8
																																					R1b1a1b1a1a2c1a4b2a1d	DC25
																																			R1b1a1b1a1a2c1a4b2b	A16
																																				R1b1a1b1a1a2c1a4b2b1	A14
																																					R1b1a1b1a1a2c1a4b2b1a	L643
																																			R1b1a1b1a1a2c1a4b2c	CTS9975
																																				R1b1a1b1a1a2c1a4b2c1	Z2186
																																					R1b1a1b1a1a2c1a4b2c1a	CTS1202.1
																																						R1b1a1b1a1a2c1a4b2c1a1	CTS9881
																																							R1b1a1b1a1a2c1a4b2c1a1a	BY412
																																							R1b1a1b1a1a2c1a4b2c1a1b	FGC32677
																																								R1b1a1b1a1a2c1a4b2c1a1b1	Y6067
																																									R1b1a1b1a1a2c1a4b2c1a1b1a	Z18128
																																									R1b1a1b1a1a2c1a4b2c1a1b1b	CTS4296
																																										R1b1a1b1a1a2c1a4b2c1a1b1b1	A2072.2
																																											R1b1a1b1a1a2c1a4b2c1a1b1b1a	Z33816.2
																																					R1b1a1b1a1a2c1a4b2c1b	A2298
																																						R1b1a1b1a1a2c1a4b2c1b1	BY2848
																																							R1b1a1b1a1a2c1a4b2c1b1a	BY2760
																																				R1b1a1b1a1a2c1a4b2c2	S891
																																				R1b1a1b1a1a2c1a4b2c3	A287
																																					R1b1a1b1a1a2c1a4b2c3a	FGC36536
																																				R1b1a1b1a1a2c1a4b2c4	FGC17551
																																			R1b1a1b1a1a2c1a4b2d	S939
																																				R1b1a1b1a1a2c1a4b2d1	DF73
																																				R1b1a1b1a1a2c1a4b2d2	BY279
																																			R1b1a1b1a1a2c1a4b3a	FGC3222
																																				R1b1a1b1a1a2c1a4b3a1	PF825.2
																																					R1b1a1b1a1a2c1a4b3a1a	FGC3221
																																						R1b1a1b1a1a2c1a4b3a1a1	FGC3251
																																		R1b1a1b1a1a2c1a4b4	S844
																																			R1b1a1b1a1a2c1a4b4a	S856
																																				R1b1a1b1a1a2c1a4b4a1	S845
																																					R1b1a1b1a1a2c1a4b4a1a	S846
																																						R1b1a1b1a1a2c1a4b4a1a1	FGC20560
																																							R1b1a1b1a1a2c1a4b4a1a1a	FGC20561
																																				R1b1a1b1a1a2c1a4b4a2	L1308
																																		R1b1a1b1a1a2c1a4b5	A494
																																			R1b1a1b1a1a2c1a4b5a	S7898
																																				R1b1a1b1a1a2c1a4b5a1	Z18132
																																					R1b1a1b1a1a2c1a4b5a1a	BY157
																																						R1b1a1b1a1a2c1a4b5a1a1	BY357
																																				R1b1a1b1a1a2c1a4b5a2	S15280.1
																																					R1b1a1b1a1a2c1a4b5a2a	FGC17449
																																		R1b1a1b1a1a2c1a4b6	BY312
																																			R1b1a1b1a1a2c1a4b6a	BY325
																																	R1b1a1b1a1a2c1a4c	CTS3386
																																		R1b1a1b1a1a2c1a4c1	S19268
																																	R1b1a1b1a1a2c1a4d	MC14
																																		R1b1a1b1a1a2c1a4d1	A2070
																																		R1b1a1b1a1a2c1a4d2	A5384
																																R1b1a1b1a1a2c1a5	DF21
																																	R1b1a1b1a1a2c1a5a	FGC3213
																																		R1b1a1b1a1a2c1a5a1	Z16532
																																			R1b1a1b1a1a2c1a5a1a	P314.2
																																				R1b1a1b1a1a2c1a5a1a1	L294.2
																																				R1b1a1b1a1a2c1a5a1a2	Z16533
																																		R1b1a1b1a1a2c1a5a2	ZZ1
																																			R1b1a1b1a1a2c1a5a2a	S3058
																																				R1b1a1b1a1a2c1a5a2a1	S424
																																					R1b1a1b1a1a2c1a5a2a1a	S426
																																						R1b1a1b1a1a2c1a5a2a1a1	CTS2187
																																							R1b1a1b1a1a2c1a5a2a1a1a	S425
																																						R1b1a1b1a1a2c1a5a2a1a2	Z36747
																																						R1b1a1b1a1a2c1a5a2a1a3	FGC3197
																																						R1b1a1b1a1a2c1a5a2a1a4	ZZ23_1. ZZ23_2
																																				R1b1a1b1a1a2c1a5a2a2	BY2719
																																					R1b1a1b1a1a2c1a5a2a2a	BY3348
																																						R1b1a1b1a1a2c1a5a2a2a1	BY2198
																																				R1b1a1b1a1a2c1a5a2a3	Y22972
																																			R1b1a1b1a1a2c1a5a2b	S5459
																																				R1b1a1b1a1a2c1a5a2b1	S5456
																																					R1b1a1b1a1a2c1a5a2b1a	Z29558
																																						R1b1a1b1a1a2c1a5a2b1a1	BY3736
																																					R1b1a1b1a1a2c1a5a2b1b	BY3300
																																	R1b1a1b1a1a2c1a5b	Z30233
																																		R1b1a1b1a1a2c1a5b1	FGC3903
																																			R1b1a1b1a1a2c1a5b1a	S280
																																				R1b1a1b1a1a2c1a5b1a1	DF25
																																					R1b1a1b1a1a2c1a5b1a1a	DF5
																																						R1b1a1b1a1a2c1a5b1a1a1	RC6
																																							R1b1a1b1a1a2c1a5b1a1a1a	FGC5780
																																								R1b1a1b1a1a2c1a5b1a1a1a1	RC5
																																									R1b1a1b1a1a2c1a5b1a1a1a1a	L658
																																										R1b1a1b1a1a2c1a5b1a1a1a1a1	FGC5787
																																											R1b1a1b1a1a2c1a5b1a1a1a1a1a	FGC5766
																																												R1b1a1b1a1a2c1a5b1a1a1a1a1a1	FGC5760
																																											R1b1a1b1a1a2c1a5b1a1a1a1a1b	F15489
																																									R1b1a1b1a1a2c1a5b1a1a1a1b	BY9595
																																										R1b1a1b1a1a2c1a5b1a1a1a1b1	BY11494
																																						R1b1a1b1a1a2c1a5b1a1a2	CTS3655
																																							R1b1a1b1a1a2c1a5b1a1a2a	L627
																																								R1b1a1b1a1a2c1a5b1a1a2a1	BY3364
																																									R1b1a1b1a1a2c1a5b1a1a2a1a	BY12309
																																							R1b1a1b1a1a2c1a5b1a1a2b	Z16539
																																								R1b1a1b1a1a2c1a5b1a1a2b1	Z16540
																																									R1b1a1b1a1a2c1a5b1a1a2b1a	BY3363
																																										R1b1a1b1a1a2c1a5b1a1a2b1a1	V75.2
																																						R1b1a1b1a1a2c1a5b1a1a3	L1403
																																							R1b1a1b1a1a2c1a5b1a1a3a	L1402
																																								R1b1a1b1a1a2c1a5b1a1a3a1	A421
																																								R1b1a1b1a1a2c1a5b1a1a3a2	A426
																																								R1b1a1b1a1a2c1a5b1a1a3a3	A459
																																							R1b1a1b1a1a2c1a5b1a1a3b	A5409
																																								R1b1a1b1a1a2c1a5b1a1a3b1	BY3455
																																									R1b1a1b1a1a2c1a5b1a1a3b1a	A418
																																										R1b1a1b1a1a2c1a5b1a1a3b1a1	BY4132
																																						R1b1a1b1a1a2c1a5b1a1a4	FGC15498
																																							R1b1a1b1a1a2c1a5b1a1a4a	L1446
																																								R1b1a1b1a1a2c1a5b1a1a4a1	L1448
																																								R1b1a1b1a1a2c1a5b1a1a4a2	FGC15499
																																					R1b1a1b1a1a2c1a5b1a1b	Z29528
																																					R1b1a1b1a1a2c1a5b1a1c	S6189
																																						R1b1a1b1a1a2c1a5b1a1c1	S6168
																																							R1b1a1b1a1a2c1a5b1a1c1a	S6185
																																			R1b1a1b1a1a2c1a5b1b	FGC9749
																																	R1b1a1b1a1a2c1a5c	S5488
																																		R1b1a1b1a1a2c1a5c1	S7200
																																			R1b1a1b1a1a2c1a5c1a	L720
																																				R1b1a1b1a1a2c1a5c1a1	BY11111
																																					R1b1a1b1a1a2c1a5c1a1a	A14025
																																					R1b1a1b1a1a2c1a5c1a1b	BY30496
																																				R1b1a1b1a1a2c1a5c1a2	FGC34092
																																			R1b1a1b1a1a2c1a5c1b	S6003
																																				R1b1a1b1a1a2c1a5c1b1	Z16541
																																					R1b1a1b1a1a2c1a5c1b1a	S6000
																																						R1b1a1b1a1a2c1a5c1b1a1	Z17557
																																							R1b1a1b1a1a2c1a5c1b1a1a	BY20276
																																						R1b1a1b1a1a2c1a5c1b1a2	Y23529
																																						R1b1a1b1a1a2c1a5c1b1a3	FGC43483
																																							R1b1a1b1a1a2c1a5c1b1a3a	FGC43477
																																								R1b1a1b1a1a2c1a5c1b1a3a1	A11730
																																									R1b1a1b1a1a2c1a5c1b1a3a1a	A11733
																																					R1b1a1b1a1a2c1a5c1b1b	BY634
																																						R1b1a1b1a1a2c1a5c1b1b1	FGC33971
																																		R1b1a1b1a1a2c1a5c2	FGC33060
																																			R1b1a1b1a1a2c1a5c2b	FGC33039
																																				R1b1a1b1a1a2c1a5c2b1	FGC33041
																																		R1b1a1b1a1a2c1a5c3	Z16294
																																			R1b1a1b1a1a2c1a5c3a	L130
																																			R1b1a1b1a1a2c1a5c3b	Z16292
																																				R1b1a1b1a1a2c1a5c3b1	Z16282
																																					R1b1a1b1a1a2c1a5c3b1a	Z16284
																																						R1b1a1b1a1a2c1a5c3b1a1	Z16289
																																							R1b1a1b1a1a2c1a5c3b1a1a	BY20009
																																								R1b1a1b1a1a2c1a5c3b1a1a1	BY20010
																																						R1b1a1b1a1a2c1a5c3b1a2	BY3301
																																							R1b1a1b1a1a2c1a5c3b1a2a	FGC63585
																																							R1b1a1b1a1a2c1a5c3b1a2b	BY30490
																																						R1b1a1b1a1a2c1a5c3b1a3	BY4011
																																							R1b1a1b1a1a2c1a5c3b1a3a	BY4005
																																								R1b1a1b1a1a2c1a5c3b1a3a1	FGC61810
																																					R1b1a1b1a1a2c1a5c3b1b	BY3302
																																						R1b1a1b1a1a2c1a5c3b1b1	BY3303
																																						R1b1a1b1a1a2c1a5c3b1b2	BY23092
																																							R1b1a1b1a1a2c1a5c3b1b2a	BY31354
																																				R1b1a1b1a1a2c1a5c3b2	A14023
																																			R1b1a1b1a1a2c1a5c3c	BY12084
																																				R1b1a1b1a1a2c1a5c3c1	FGC58881
																																			R1b1a1b1a1a2c1a5c3d	S19938
																																				R1b1a1b1a1a2c1a5c3d1	BY35673
																																			R1b1a1b1a1a2c1a5c3e	FGC42958
																																			R1b1a1b1a1a2c1a5c3f	BY32772
																																			R1b1a1b1a1a2c1a5c3g	A16411
																																		R1b1a1b1a1a2c1a5c4	FGC11358
																																			R1b1a1b1a1a2c1a5c4a	L1336
																																				R1b1a1b1a1a2c1a5c4a1	FGC11336
																																				R1b1a1b1a1a2c1a5c4a2	A933
																																					R1b1a1b1a1a2c1a5c4a2a	A935
																																						R1b1a1b1a1a2c1a5c4a2a1	A932
																																				R1b1a1b1a1a2c1a5c4a3	BY15552
																																			R1b1a1b1a1a2c1a5c4b	BY19855
																																		R1b1a1b1a1a2c1a5c5	A6484
																																			R1b1a1b1a1a2c1a5c5a	BY8982
																																				R1b1a1b1a1a2c1a5c5a1	BY8980
																																		R1b1a1b1a1a2c1a5c6	BY12129
																																	R1b1a1b1a1a2c1a5d	S971
																																		R1b1a1b1a1a2c1a5d1	FGC23381
																																			R1b1a1b1a1a2c1a5d1a	FGC23375
																																		R1b1a1b1a1a2c1a5d2	FGC19706
																																		R1b1a1b1a1a2c1a5d3	A51
																																	R1b1a1b1a1a2c1a5e	BY9405
																																R1b1a1b1a1a2c1a6	FGC5494
																																	R1b1a1b1a1a2c1a6a	FGC5496
																																		R1b1a1b1a1a2c1a6a1	S1088
																																			R1b1a1b1a1a2c1a6a1a	CTS2457.2
																															R1b1a1b1a1a2c1b	CTS300
																																R1b1a1b1a1a2c1b1	CTS6919
																																	R1b1a1b1a1a2c1b1a	A92
																																R1b1a1b1a1a2c1b2	BY592
																																	R1b1a1b1a1a2c1b2a	Z16245
																																		R1b1a1b1a1a2c1b2a1	Z16246
																																			R1b1a1b1a1a2c1b2a1a	BY3070
																																				R1b1a1b1a1a2c1b2a1a1	BY3004
																																					R1b1a1b1a1a2c1b2a1a1a	FGC36421
																																			R1b1a1b1a1a2c1b2a1b	BY2583
																																	R1b1a1b1a1a2c1b2b	A7811
																															R1b1a1b1a1a2c1c	A5846
																																R1b1a1b1a1a2c1c1	A5840
																															R1b1a1b1a1a2c1d	A7905
																													R1b1a1b1a1a2d	L238
																														R1b1a1b1a1a2d1	Z2245
																															R1b1a1b1a1a2d1a	Z2247
																																R1b1a1b1a1a2d1a1	CTS11638
																													R1b1a1b1a1a2e	DF19
																														R1b1a1b1a1a2e1	DF88
																															R1b1a1b1a1a2e1a	L644
																															R1b1a1b1a1a2e1b	L1199
																															R1b1a1b1a1a2e1c	L719
																														R1b1a1b1a1a2e2	S233
																													R1b1a1b1a1a2f	DF99
																													R1b1a1b1a1a2g	Y18209
																													R1b1a1b1a1a2h	A9063
																												R1b1a1b1a1a3	AM01876
																													R1b1a1b1a1a3a	AM01877
																														R1b1a1b1a1a3a1	S14328
																														R1b1a1b1a1a3a2	S11481
																													R1b1a1b1a1a3b	A8039
																														R1b1a1b1a1a3b1	A8040
																															R1b1a1b1a1a3b1a	A8041
																															R1b1a1b1a1a3b1b	A8045
																												R1b1a1b1a1a4	A8051
																													R1b1a1b1a1a4a	A8055
																													R1b1a1b1a1a4b	FGC37083
																										R1b1a1b1a2	PF7589
																											R1b1a1b1a2a	CTS6889
																												R1b1a1b1a2a1	CTS11824
																									R1b1a1b1b	CTS1078
																										R1b1a1b1b1	L584
																											R1b1a1b1b1a	L943
																										R1b1a1b1b2	L277.1
																										R1b1a1b1b3	Z2106
																											R1b1a1b1b3a	CTS1843
																												R1b1a1b1b3a1	CTS7822
																													R1b1a1b1b3a1a	CTS7556
																														R1b1a1b1b3a1a1	CTS9219
																															R1b1a1b1b3a1a1a	Y5587
																															R1b1a1b1b3a1a1b	BY250
																															R1b1a1b1b3a1a1c	A1777
																												R1b1a1b1b3a2	Y14416
																								R1b1a1b2	PF7558
																									R1b1a1b2a	GG480
																										R1b1a1b2a1	Y31335
																											R1b1a1b2a1a	BY16680
																											R1b1a1b2a1b	A11710
																												R1b1a1b2a1b1	A15807
																													R1b1a1b2a1b1a	A15808
																													R1b1a1b2a1b1b	A15809
																														R1b1a1b2a1b1b1	BY16696
																										R1b1a1b2a2	Z29758
																											R1b1a1b2a2a	PF7566
																											R1b1a1b2a2b	GG700
																												R1b1a1b2a2b1	V3286
																												R1b1a1b2a2b2	GG445
																												R1b1a1b2a2b3	FGC42003
																													R1b1a1b2a2b3a	Y31465
																												R1b1a1b2a2b4	PH2558
																													R1b1a1b2a2b4a	PH4238
																									R1b1a1b2b	FGC31923
																						R1b1a2	V1636
																							R1b1a2a	V1274
																								R1b1a2a1	BY15381
																									R1b1a2a1a	CTS5330
																									R1b1a2a1c	BY46176
																							R1b1a2b	Y83069
																								R1b1a2b1	BY139812
																					R1b1b	PF6279
																						R1b1b1	M18
																						R1b1b2	FGC21014
																							R1b1b2a	FGC21034
																								R1b1b2a1	V35
																								R1b1b2a2	FGC20970
																									R1b1b2a2a	SK2071
																										R1b1b2a2a1	V69
																											R1b1b2a2a1a	FGC39691
																												R1b1b2a2a1a1	Y24712
																										R1b1b2a2a2	FGC42146
																									R1b1b2a2b	FGC20973
																										R1b1b2a2b1	FGC20980
																											R1b1b2a2b1a	FGC21047
																												R1b1b2a2b1a1	FGC21060
																				R1b2	PH155
																					R1b2a	M335
																					R1b2b	PH200
																						R1b2b1	PH1578
																							R1b2b1a	BY14576
																								R1b2b1a1	Y32791
																						R1b2b2	Y92825
																		R2	M479
																			R2a	M124
																				R2a1	L263
																				R2a2	P267
																					R2a2a	FGC13203
																						R2a2a1	FGC13188
																							R2a2a1a	F1092
																								R2a2a1a1	F1758
																									R2a2a1a1a	FGC13184
																										R2a2a1a1a1	FGC13201
																											R2a2a1a1a1a	BY14642
																												R2a2a1a1a1a1	FGC13235
																					R2a2b	Y12100
																						R2a2b1	Y8763
																							R2a2b1a	L1069
																							R2a2b1b	FGC17608
																								R2a2b1b1	FGC17611
																									R2a2b1b1a	M2091
																										R2a2b1b1a1	FGC17661
																										R2a2b1b1a2	BY21669
																								R2a2b1b2	FGC18148
																									R2a2b1b2a	SK2142
																										R2a2b1b2a1	Y1377
																											R2a2b1b2a1a	Y1379
																												R2a2b1b2a1a1	Y1383
																									R2a2b1b2b	FGC18147
																										R2a2b1b2b1	L723
																											R2a2b1b2b1a	L725
																												R2a2b1b2b1a1	L724
																										R2a2b1b2b2	Y1283
																											R2a2b1b2b2a	SK2155
																												R2a2b1b2b2a1	Y1331
																													R2a2b1b2b2a1a	Y20020
																													R2a2b1b2b2a1b	Y26635
																														R2a2b1b2b2a1b1	Y34210
																												R2a2b1b2b2a2	Y28598
																											R2a2b1b2b2b	Z29162
																												R2a2b1b2b2b1	V4082
																													R2a2b1b2b2b1a	YP5340
																													R2a2b1b2b2b1b	Y28589
																												R2a2b1b2b2b2	Y1292
																													R2a2b1b2b2b2a	Y1296
																														R2a2b1b2b2b2a1	Y1303
																															R2a2b1b2b2b2a1a	Y1325
																																R2a2b1b2b2b2a1a1	Y2174
																										R2a2b1b2b3	V3467
																											R2a2b1b2b3a	Y1357
																											R2a2b1b2b3b	V2434
																												R2a2b1b2b3b1	Y1371
																												R2a2b1b2b3b2	L294.1
																													R2a2b1b2b3b2a	V1024
																														R2a2b1b2b3b2a1	FGC18152
																															R2a2b1b2b3b2a1a	Y17474
																															R2a2b1b2b3b2a1b	FGC7156
																									R2a2b1b2c	Y26630
																										R2a2b1b2c1	Y28641
																						R2a2b2	FGC46676
																							R2a2b2a	FGC51802
																							R2a2b2b	BY14309
																								R2a2b2b1	L288.1
																			R2b	FGC21706
																				R2b1	FGC50339
																					R2b1a	FGC50273
																						R2b1a1	FGC50249
"""


class GenotypeStatus:
    """Genotype status enumeration"""
    DERIVED = "derived"      # Derived
    ANCESTRAL = "ancestral"  # Ancestral
    HETEROZYGOUS = "het"     # Heterozygous
    MISSING = "missing"      # Missing
    MISMATCH = "mismatch"    # Mismatch


class HaplogroupClassifier:
    """
    Y-chromosome Haplogroup Classifier
    
    Core logic: simulates manual classification
    1. Traverse phylogenetic tree from root downward
    2. Check defining SNP status at each node
    3. Derived → enter, Ancestral → stop, Missing → check downstream
    4. Validate ancestor chain consistency
    """
    
    # Terminal main branches (includes all major global branches)
    # Note: F is ancestor of GHIJK, most cases should be subdivided to G/H/I/J/K downstream
    # True basal F is rare, handled by special logic
    MAIN_BRANCHES = ['A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'L', 'N', 'O', 'Q', 'R', 'T']
    
    # Intermediate nodes (not used as final result unless cannot subdivide)
    INTERMEDIATE_NODES = {
        'K', 'NO', 'NO1', 'P', 'P1', 'LT', 'F', 'CF', 'DE', 'CT', 
        'BT', 'GHIJK', 'HIJK', 'IJK', 'IJ', 'A1b', 'A1b1',
        'A0-T', 'A00-T', 'A000-T'
    }
    
    # Main branch defining SNPs
    MAIN_BRANCH_DEFINING_SNPS = {
        'A': 'L1085',  # A0-TDefining SNP，Covers all A lineages
        'B': 'M60',    # BHaplogroup defining SNP
        'C': 'M130', 'D': 'M174', 'E': 'M96',
        'G': 'M201', 'H': 'L901', 'I': 'M170',
        'J': 'M304', 'L': 'M20', 'N': 'M231',
        'O': 'M175',  # indel，May not be detectable
        'Q': 'M242', 'R': 'M207', 'T': 'M184'
    }
    
    # FHaplogroupSpecial handling（basal FDetection）
    F_DEFINING_SNP = 'M89'
    
    def __init__(self, isogg_file, het_mode="moderate"):
        """
        Initialize classifier
        
        het_mode:
            strict:   Heterozygous treated as missing (most conservative)
            moderate: Heterozygous treated as derived but marked (default)
            lenient:  Heterozygous fully treated as derived
        """
        self.het_mode = het_mode
        
        print(f"[1] Parsing phylogenetic tree...")
        self.phylo_tree, self.node_snps = self._parse_tree(OFFICIAL_TREE)
        print(f"    Node count: {len(self.phylo_tree)}")
        
        print(f"\n[2] Reading ISOGG index: {isogg_file}")
        self.snp_to_info = {}      # snp_name -> (pos, ref, alt)
        self.pos_to_haplo = defaultdict(list)  # pos -> [(haplo, snp, ref, alt), ...]
        self._load_isogg(isogg_file)
        
        # Build node to SNP position mapping
        self._build_node_snp_positions()
        
        # Check main branch SNPs
        self._check_main_branch_snps()
    
    def _parse_tree(self, tree_text):
        """Parse indented format phylogenetic tree"""
        phylo_tree = {}  # child -> parent
        node_snps = {}   # node -> [snp1, snp2, ...]
        
        lines = tree_text.strip().split('\n')
        path_stack = []  # [(indent, node_name), ...]
        
        for line in lines:
            if not line.strip():
                continue
            
            # Calculate indentation
            indent = 0
            for char in line:
                if char == '\t':
                    indent += 1
                else:
                    break
            
            content = line.strip()
            parts = content.split('\t')
            
            if len(parts) < 1:
                continue
            
            node_name = parts[0].strip()
            
            # Skip root node and invalid lines
            if node_name in ['Y', 'Root', ''] or 'see ' in node_name.lower():
                continue
            
            # Extract SNPs
            snps = []
            if len(parts) >= 2:
                snp_text = parts[1].strip()
                for snp_part in snp_text.split(','):
                    snp_part = snp_part.strip()
                    if snp_part and not snp_part.startswith('('):
                        for alias in snp_part.split('/'):
                            alias = alias.strip().rstrip('~^*')
                            if alias:
                                snps.append(alias)
            
            node_snps[node_name] = snps
            
            # Maintain path stack
            while path_stack and path_stack[-1][0] >= indent:
                path_stack.pop()
            
            if path_stack:
                parent = path_stack[-1][1]
                phylo_tree[node_name] = parent
            
            path_stack.append((indent, node_name))
        
        return phylo_tree, node_snps
    
    def _load_isogg(self, isogg_file):
        """Load ISOGG index file"""
        count = 0
        
        with open(isogg_file, 'r', encoding='utf-8', errors='ignore') as f:
            header = f.readline()
            sep = '\t' if '\t' in header else ','
            
            for line in f:
                parts = line.strip().split(sep)
                if len(parts) < 7:
                    continue
                
                snp_name = parts[0].strip().rstrip('~^*')
                haplo = parts[1].strip().rstrip('~^*')
                pos_str = parts[4].strip()
                mutation = parts[6].strip() if len(parts) > 6 else ""
                
                if not pos_str or '..' in pos_str:
                    continue
                
                try:
                    pos = int(pos_str)
                except:
                    continue
                
                ref, alt = '', ''
                if '->' in mutation:
                    try:
                        ref, alt = mutation.split('->')[:2]
                        ref, alt = ref.strip().upper(), alt.strip().upper()
                    except:
                        pass
                
                self.snp_to_info[snp_name] = (pos, ref, alt)
                self.pos_to_haplo[pos].append((haplo, snp_name, ref, alt))
                count += 1
                
                # Process aliases
                if len(parts) > 2:
                    for alt_name in parts[2].split(';'):
                        alt_name = alt_name.strip().rstrip('~^*')
                        if alt_name and alt_name not in self.snp_to_info:
                            self.snp_to_info[alt_name] = (pos, ref, alt)
        
        print(f"    Loaded SNPs: {count}")
    
    def _build_node_snp_positions(self):
        """Build SNP position mapping for tree nodes"""
        self.node_snp_positions = {}
        
        for node, snps in self.node_snps.items():
            positions = []
            for snp in snps:
                if snp in self.snp_to_info:
                    pos, ref, alt = self.snp_to_info[snp]
                    positions.append((snp, pos, ref, alt))
            self.node_snp_positions[node] = positions
        
        nodes_with_pos = sum(1 for p in self.node_snp_positions.values() if p)
        print(f"    {nodes_with_pos}/{len(self.node_snps)} tree nodes have SNP positions")
    
    def _check_main_branch_snps(self):
        """Check availability of main branch defining SNPs"""
        print(f"\n    Main branch SNP check:")
        
        for branch in self.MAIN_BRANCHES:
            snp = self.MAIN_BRANCH_DEFINING_SNPS.get(branch, '')
            if snp and snp in self.snp_to_info:
                pos = self.snp_to_info[snp][0]
                print(f"      {branch}: ✓ {snp} @ {pos}")
            else:
                print(f"      {branch}: ✗ {snp} (Needs inference from downstream)")
    
    def get_ancestors(self, node):
        """Get all ancestors of node (nearest to farthest)"""
        ancestors = []
        current = node
        seen = set()
        
        while current in self.phylo_tree and current not in seen:
            seen.add(current)
            parent = self.phylo_tree[current]
            ancestors.append(parent)
            current = parent
        
        return ancestors
    
    def get_depth(self, node):
        """Get node depth in tree"""
        return len(self.get_ancestors(node))
    
    def get_main_branch(self, haplo):
        """Extract main branch from haplogroup name"""
        if not haplo:
            return None
        
        # Direct match
        for branch in self.MAIN_BRANCHES:
            if haplo == branch:
                return branch
            if haplo.startswith(branch) and len(haplo) > 1:
                next_char = haplo[len(branch)] if len(haplo) > len(branch) else ''
                if next_char.isdigit() or next_char == '':
                    return branch
        
        # Check if intermediate node
        if haplo in self.INTERMEDIATE_NODES:
            return None
        
        # First letter match
        first = haplo[0] if haplo else None
        if first in self.MAIN_BRANCHES:
            return first
        
        return None
    
    def check_genotype(self, gt, vcf_ref, vcf_alt, exp_ref, exp_alt):
        """
        Check genotype status
        
        Returns: GenotypeStatus
        """
        if gt in ['./.', '.', '']:
            return GenotypeStatus.MISSING
        
        vcf_ref = vcf_ref.upper()
        vcf_alt = vcf_alt.upper() if vcf_alt != '.' else ''
        
        is_het = gt in ['0/1', '0|1', '1/0', '1|0']
        is_hom_alt = gt in ['1/1', '1|1', '1']
        is_hom_ref = gt in ['0/0', '0|0', '0']
        
        # Has explicit ref/alt information
        if exp_ref and exp_alt:
            # Forward match: VCF REF=ancestral, ALT=derived
            if vcf_ref == exp_ref.upper() and vcf_alt == exp_alt.upper():
                if is_hom_alt:
                    return GenotypeStatus.DERIVED
                elif is_het:
                    return GenotypeStatus.HETEROZYGOUS
                elif is_hom_ref:
                    return GenotypeStatus.ANCESTRAL
            # Reverse match: VCF REF=derived, ALT=ancestral
            elif vcf_ref == exp_alt.upper() and vcf_alt == exp_ref.upper():
                if is_hom_ref:
                    return GenotypeStatus.DERIVED
                elif is_het:
                    return GenotypeStatus.HETEROZYGOUS
                elif is_hom_alt:
                    return GenotypeStatus.ANCESTRAL
        
        # When no explicit info, infer from VCF
        if is_hom_alt:
            return GenotypeStatus.DERIVED
        elif is_het:
            return GenotypeStatus.HETEROZYGOUS
        elif is_hom_ref:
            return GenotypeStatus.ANCESTRAL
        
        return GenotypeStatus.MISMATCH
    
    def is_derived(self, status):
        """
        Determine if treated as derived based on het_mode
        
        Returns: True/False/None (NoneIndicates missing)
        """
        if status == GenotypeStatus.DERIVED:
            return True
        elif status == GenotypeStatus.ANCESTRAL:
            return False
        elif status == GenotypeStatus.HETEROZYGOUS:
            if self.het_mode == "strict":
                return None
            else:  # moderate or lenient
                return True
        else:  # MISSING or MISMATCH
            return None
    
    def check_node_status(self, node, sample_geno):
        """
        Check status of single node
        
        Returns: {
            'status': 'derived'/'ancestral'/'missing'/'mixed',
            'derived_snps': [...],
            'ancestral_snps': [...],
            'het_snps': [...],
            'missing_snps': [...]
        }
        """
        positions = self.node_snp_positions.get(node, [])
        
        derived_snps = []
        ancestral_snps = []
        het_snps = []
        missing_snps = []
        
        for snp, pos, ref, alt in positions:
            if pos not in sample_geno:
                missing_snps.append(snp)
                continue
            
            gt, vcf_ref, vcf_alt = sample_geno[pos]
            status = self.check_genotype(gt, vcf_ref, vcf_alt, ref, alt)
            
            if status == GenotypeStatus.DERIVED:
                derived_snps.append(snp)
            elif status == GenotypeStatus.ANCESTRAL:
                ancestral_snps.append(snp)
            elif status == GenotypeStatus.HETEROZYGOUS:
                het_snps.append(snp)
            else:
                missing_snps.append(snp)
        
        # Overall determination
        if derived_snps:
            overall = 'derived'
        elif ancestral_snps and not derived_snps:
            overall = 'ancestral'
        elif het_snps and not derived_snps and not ancestral_snps:
            overall = 'het'
        else:
            overall = 'missing'
        
        return {
            'status': overall,
            'derived_snps': derived_snps,
            'ancestral_snps': ancestral_snps,
            'het_snps': het_snps,
            'missing_snps': missing_snps
        }
    
    def classify_sample(self, sample, sample_geno):
        """
        Classify single sample
        
        Core logic：
        1. Collect all derived haplogroups
        2. Group by main branch
        3. Validate ancestor chain for each candidate
        4. Select deepest validated node
        """
        het_count = 0
        
        # ========================================
        # Step 1: Collect all derived SNPs/haplogroups
        # ========================================
        derived_haplos = defaultdict(list)  # haplo -> [(snp, pos), ...]
        all_node_status = {}  # Record status of all checked nodes
        
        for pos, entries in self.pos_to_haplo.items():
            if pos not in sample_geno:
                continue
            
            gt, vcf_ref, vcf_alt = sample_geno[pos]
            
            # Count heterozygous
            if gt in ['0/1', '0|1', '1/0', '1|0']:
                het_count += 1
            
            for haplo, snp_name, ref, alt in entries:
                status = self.check_genotype(gt, vcf_ref, vcf_alt, ref, alt)
                
                # Record node status
                if haplo not in all_node_status:
                    all_node_status[haplo] = {
                        'derived': [], 'ancestral': [], 'het': [], 'missing': []
                    }
                
                if status == GenotypeStatus.DERIVED:
                    all_node_status[haplo]['derived'].append(snp_name)
                elif status == GenotypeStatus.ANCESTRAL:
                    all_node_status[haplo]['ancestral'].append(snp_name)
                elif status == GenotypeStatus.HETEROZYGOUS:
                    all_node_status[haplo]['het'].append(snp_name)
                
                # Determine if derived
                if self.is_derived(status):
                    derived_haplos[haplo].append((snp_name, pos))
        
        if not derived_haplos:
            return self._make_result(
                sample, 'Unknown', 'Unknown', 0, 0, het_count, 
                [], 'No derived SNPs', {}
            )
        
        # ========================================
        # Step 2: Group by main branch
        # ========================================
        branch_candidates = defaultdict(list)  # branch -> [(haplo, n_snps, snps), ...]
        
        for haplo, snps in derived_haplos.items():
            branch = self.get_main_branch(haplo)
            if branch:
                branch_candidates[branch].append({
                    'haplo': haplo,
                    'n_snps': len(snps),
                    'snps': [s[0] for s in snps],
                    'depth': self.get_depth(haplo)
                })
        
        if not branch_candidates:
            return self._make_result(
                sample, 'Unknown', 'Unknown', 0, 0, het_count,
                [], 'Cannot determine main branch', {}
            )
        
        # ========================================
        # Step 3: Determine main branch
        # ========================================
        # Check main branch defining SNPs
        branch_direct_status = {}
        for branch in self.MAIN_BRANCHES:
            snp = self.MAIN_BRANCH_DEFINING_SNPS.get(branch)
            if snp and snp in self.snp_to_info:
                pos, ref, alt = self.snp_to_info[snp]
                if pos in sample_geno:
                    gt, vcf_ref, vcf_alt = sample_geno[pos]
                    status = self.check_genotype(gt, vcf_ref, vcf_alt, ref, alt)
                    branch_direct_status[branch] = status
                else:
                    branch_direct_status[branch] = GenotypeStatus.MISSING
            else:
                branch_direct_status[branch] = GenotypeStatus.MISSING
        
        # Select main branch: prioritize defining SNPs, then downstream evidence
        main_branch = None
        main_branch_note = ""
        
        # First find those with derived defining SNPs
        for branch in self.MAIN_BRANCHES:
            status = branch_direct_status.get(branch)
            if self.is_derived(status) and branch in branch_candidates:
                main_branch = branch
                main_branch_note = f"{self.MAIN_BRANCH_DEFINING_SNPS[branch]} derived"
                break
        
        # Not found, infer from downstream evidence
        if not main_branch:
            # Exclude branches with clear ancestral defining SNPs
            excluded = set()
            for branch, status in branch_direct_status.items():
                if status == GenotypeStatus.ANCESTRAL:
                    excluded.add(branch)
            
            # Select one with most downstream evidence
            best_branch = None
            best_count = 0
            
            for branch, candidates in branch_candidates.items():
                if branch in excluded:
                    continue
                total_snps = sum(c['n_snps'] for c in candidates)
                if total_snps > best_count:
                    best_count = total_snps
                    best_branch = branch
            
            if best_branch:
                main_branch = best_branch
                main_branch_note = f"inferred from {best_count} downstream SNPs"
        
        # Special handling: detect basal F
        # If no main branch found but M89 derived, may be rare basal F
        if not main_branch:
            if self.F_DEFINING_SNP in self.snp_to_info:
                f_pos = self.snp_to_info[self.F_DEFINING_SNP][0]
                if f_pos in sample_geno:
                    gt, vcf_ref, vcf_alt = sample_geno[f_pos]
                    f_ref, f_alt = self.snp_to_info[self.F_DEFINING_SNP][1], self.snp_to_info[self.F_DEFINING_SNP][2]
                    f_status = self.check_genotype(gt, vcf_ref, vcf_alt, f_ref, f_alt)
                    if self.is_derived(f_status):
                        # M89Derived but no downstream branch evidence = basal F
                        main_branch = 'F'
                        main_branch_note = "basal F (M89 derived, no downstream)"
                        # Add F to branch_candidates for subsequent processing
                        if 'F' not in branch_candidates:
                            branch_candidates['F'] = [{'haplo': 'F', 'n_snps': 1, 'snps': ['M89'], 'depth': 0}]
        
        if not main_branch:
            return self._make_result(
                sample, 'Unknown', 'Unknown', 0, 0, het_count,
                [], 'Cannot determine main branch（Insufficient or conflicting evidence）', {}
            )
        
        # ========================================
        # Step 4: Subdivide within main branch with ancestor chain validation
        # ========================================
        candidates = branch_candidates[main_branch]
        
        # Sort by depth and SNP count
        candidates.sort(key=lambda x: (x['depth'], x['n_snps']), reverse=True)
        
        # Validate ancestor chain for each candidate
        best_haplo = main_branch
        best_n_snps = 0
        best_evidence = []
        ancestor_conflicts = []
        
        for cand in candidates:
            haplo = cand['haplo']
            ancestors = self.get_ancestors(haplo)
            
            # Only check ancestors within main branch
            relevant_ancestors = []
            for anc in ancestors:
                if anc == main_branch:
                    break
                if self.get_main_branch(anc) == main_branch or anc == main_branch:
                    relevant_ancestors.append(anc)
            
            # Check ancestor chain
            chain_valid = True
            conflicts = []
            
            for anc in relevant_ancestors:
                if anc in all_node_status:
                    anc_info = all_node_status[anc]
                    has_derived = len(anc_info['derived']) > 0
                    has_ancestral = len(anc_info['ancestral']) > 0
                    has_het = len(anc_info['het']) > 0
                    
                    # Key: if ancestor is clearly ancestral with no derived evidence
                    if has_ancestral and not has_derived:
                        # Can heterozygous save it?
                        if has_het and self.het_mode in ['moderate', 'lenient']:
                            pass  # Acceptable but marked
                        else:
                            chain_valid = False
                            conflicts.append(f"{anc}(ancestral)")
            
            if chain_valid:
                best_haplo = haplo
                best_n_snps = cand['n_snps']
                best_evidence = cand['snps'][:5]
                break  # Found first deepest validated node
            else:
                ancestor_conflicts.extend(conflicts)
        
        # If all candidates have conflicts，Select one with fewest conflicts
        if best_haplo == main_branch and candidates:
            # Re-evaluate, select one with fewest conflicts
            min_conflicts = float('inf')
            for cand in candidates:
                haplo = cand['haplo']
                ancestors = self.get_ancestors(haplo)
                
                conflict_count = 0
                for anc in ancestors:
                    if anc in all_node_status:
                        anc_info = all_node_status[anc]
                        if anc_info['ancestral'] and not anc_info['derived']:
                            conflict_count += 1
                
                if conflict_count < min_conflicts:
                    min_conflicts = conflict_count
                    best_haplo = haplo
                    best_n_snps = cand['n_snps']
                    best_evidence = cand['snps'][:5]
        
        # ========================================
        # Step 5: Calculate confidence
        # ========================================
        direct_status = branch_direct_status.get(main_branch, GenotypeStatus.MISSING)
        
        # Base confidence
        if self.is_derived(direct_status):
            # Main branch defining SNP directly derived
            if best_n_snps >= 5:
                confidence = 0.95
            elif best_n_snps >= 3:
                confidence = 0.90
            else:
                confidence = 0.85
        else:
            # Inferred from downstream
            if best_n_snps >= 10:
                confidence = 0.85
            elif best_n_snps >= 5:
                confidence = 0.75
            elif best_n_snps >= 3:
                confidence = 0.65
            else:
                confidence = 0.50
        
        # Penalty factors
        if ancestor_conflicts:
            confidence -= 0.15
        if het_count > 10:
            confidence -= 0.05
        
        confidence = max(0.0, min(1.0, confidence))
        
        # Build diagnostic information
        diagnostics = {
            'main_branch_status': str(direct_status),
            'main_branch_note': main_branch_note,
            'ancestor_conflicts': ancestor_conflicts[:3] if ancestor_conflicts else [],
            'total_derived_haplos': len(derived_haplos),
            'candidates_in_branch': len(candidates)
        }
        
        return self._make_result(
            sample, main_branch, best_haplo, best_n_snps, 
            confidence, het_count, best_evidence, '', diagnostics
        )
    
    def _make_result(self, sample, main_branch, haplo, n_snps, 
                     confidence, het_count, evidence, note, diagnostics):
        """Build result dictionary"""
        return {
            'sample': sample,
            'main_branch': main_branch,
            'haplogroup': haplo,
            'n_snps': n_snps,
            'confidence': confidence,
            'het_count': het_count,
            'evidence': evidence,
            'note': note,
            'diagnostics': diagnostics
        }


def load_vcf(vcf_file):
    """Load VCF file"""
    opener = open
    mode = 'r'
    if vcf_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    
    samples = []
    sample_geno = defaultdict(dict)
    y_count = 0
    
    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = parts[9:]
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            
            chrom = parts[0]
            if chrom not in ['Y', 'chrY', '24', 'chr24']:
                continue
            
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            for i, sample in enumerate(samples):
                gt = parts[9 + i].split(':')[0]
                sample_geno[sample][pos] = (gt, ref, alt)
            
            y_count += 1
    
    return samples, sample_geno, y_count


def print_summary(results, het_mode):
    """Print summary statistics"""
    print("\n" + "=" * 70)
    print("Classification Results")
    print("=" * 70)
    
    print(f"\nHeterozygosity handling modes: {het_mode}")
    
    # Confidence distribution
    high = sum(1 for r in results if r['confidence'] >= 0.8)
    med = sum(1 for r in results if 0.5 <= r['confidence'] < 0.8)
    low = sum(1 for r in results if r['confidence'] < 0.5)
    
    print(f"\n【Confidence】")
    print(f"   High(>=0.8): {high} ({high/len(results)*100:.1f}%)")
    print(f"   Medium(0.5-0.8): {med} ({med/len(results)*100:.1f}%)")
    print(f"   Low(<0.5): {low} ({low/len(results)*100:.1f}%)")
    
    # Main branch distribution
    main_counts = defaultdict(int)
    for r in results:
        main_counts[r['main_branch']] += 1
    
    print(f"\n【Main branch】")
    for branch in sorted(main_counts.keys(), key=lambda x: -main_counts[x]):
        count = main_counts[branch]
        print(f"   {branch}: {count} ({count/len(results)*100:.1f}%)")
    
    # Detailed haplogroups
    sub_counts = defaultdict(int)
    for r in results:
        sub_counts[r['haplogroup']] += 1
    
    print(f"\n【Detailed haplogroups (top 15)】")
    for haplo in sorted(sub_counts.keys(), key=lambda x: -sub_counts[x])[:15]:
        count = sub_counts[haplo]
        print(f"   {haplo:30s}: {count:4d} ({count/len(results)*100:.1f}%)")
    
    # Ancestor chain conflict statistics
    conflict_count = sum(1 for r in results 
                         if r.get('diagnostics', {}).get('ancestor_conflicts'))
    if conflict_count > 0:
        print(f"\n【Warning】")
        print(f"   {conflict_count} samples have ancestor chain conflicts (confidence reduced)")


def save_results(results, output_dir, het_mode):
    """Save results"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Main result file
    result_file = os.path.join(output_dir, 'Y_haplogroup_result_v1.0.txt')
    with open(result_file, 'w') as f:
        f.write("Sample\tMain_Branch\tHaplogroup\tConfidence\tN_SNPs\t")
        f.write("Het_Count\tMain_Status\tEvidence\tConflicts\n")
        
        for r in results:
            evidence = ','.join(r['evidence'][:5]) if r['evidence'] else '-'
            diag = r.get('diagnostics', {})
            main_status = diag.get('main_branch_status', '-')
            conflicts = ','.join(diag.get('ancestor_conflicts', [])) or '-'
            
            f.write(f"{r['sample']}\t{r['main_branch']}\t{r['haplogroup']}\t")
            f.write(f"{r['confidence']:.3f}\t{r['n_snps']}\t{r['het_count']}\t")
            f.write(f"{main_status}\t{evidence}\t{conflicts}\n")
    
    # Summary file
    summary_file = os.path.join(output_dir, 'Y_haplogroup_summary_v1.0.txt')
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write(f"Y-chromosomeHaplogroup Classification Report v{VERSION}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Sample count: {len(results)}\n")
        f.write(f"Heterozygosity handling modes: {het_mode}\n\n")
        
        # Confidence
        high = sum(1 for r in results if r['confidence'] >= 0.8)
        med = sum(1 for r in results if 0.5 <= r['confidence'] < 0.8)
        low = sum(1 for r in results if r['confidence'] < 0.5)
        
        f.write("Confidence distribution:\n")
        f.write(f"  High(>=0.8): {high} ({high/len(results)*100:.1f}%)\n")
        f.write(f"  Medium(0.5-0.8): {med} ({med/len(results)*100:.1f}%)\n")
        f.write(f"  Low(<0.5): {low} ({low/len(results)*100:.1f}%)\n\n")
        
        # Main branch
        main_counts = defaultdict(int)
        for r in results:
            main_counts[r['main_branch']] += 1
        
        f.write("Main branch distribution:\n")
        for branch in sorted(main_counts.keys(), key=lambda x: -main_counts[x]):
            count = main_counts[branch]
            f.write(f"  {branch}: {count} ({count/len(results)*100:.1f}%)\n")
        
        # Subdivide
        sub_counts = defaultdict(int)
        for r in results:
            sub_counts[r['haplogroup']] += 1
        
        f.write("\nDetailed haplogroups:\n")
        for haplo in sorted(sub_counts.keys(), key=lambda x: -sub_counts[x])[:30]:
            count = sub_counts[haplo]
            f.write(f"  {haplo}: {count} ({count/len(results)*100:.1f}%)\n")
    
    print(f"\n[5] Results saved:")
    print(f"    {result_file}")
    print(f"    {summary_file}")
    
    return result_file, summary_file


def main():
    parser = argparse.ArgumentParser(
        description=f'Y-chromosome Haplogroup Classifier v{VERSION}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Classification logic description:
  This tool simulates manual classification, traversing phylogenetic tree from root:
  1. Check defining SNP status at each node
  2. Derived → enter node, continue downward
  3. Ancestral → stop, do not enter branch  
  4. Missing → Check downstream forDerivedEvidence（Can be inferred）
  5. Validate ancestor chain consistency, exclude conflicting paths

Heterozygosity handling modes:
  strict:   Heterozygous treated as missing (most conservative)
  moderate: Heterozygous treated as derived but marked (default)
  lenient:  Heterozygous fully treated as derived
        '''
    )
    
    parser.add_argument('-v', '--vcf', required=True, help='VCF file')
    parser.add_argument('-i', '--isogg', required=True, help='ISOGG indexdata.csv')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--het-mode', choices=['strict', 'moderate', 'lenient'],
                        default='moderate', help='Heterozygosity handling mode (Default: moderate)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"Error: {args.vcf} Does not exist")
        sys.exit(1)
    if not os.path.exists(args.isogg):
        print(f"Error: {args.isogg} Does not exist")
        sys.exit(1)
    
    print("=" * 70)
    print(f"Y-chromosome Haplogroup Classifier v{VERSION}")
    print("=" * 70)
    
    # Initialize classifier
    classifier = HaplogroupClassifier(args.isogg, args.het_mode)
    
    # LoadVCF
    print(f"\n[3] Reading VCF: {args.vcf}")
    samples, sample_geno, y_count = load_vcf(args.vcf)
    print(f"    Sample count: {len(samples)}")
    print(f"    Y-chromosome sites: {y_count}")
    
    if not samples:
        print("Error: No samples")
        sys.exit(1)
    
    # Classification
    print(f"\n[4] Classifying... (Heterozygosity mode: {args.het_mode})")
    results = []
    
    for i, sample in enumerate(samples):
        if (i + 1) % 50 == 0 or i == len(samples) - 1:
            print(f"    {i + 1}/{len(samples)}", end='\r')
        
        result = classifier.classify_sample(sample, sample_geno[sample])
        results.append(result)
    
    print("\n    Complete!♡(*´∀｀*)人(*´∀｀*)♡")
    
    # Output
    print_summary(results, args.het_mode)
    save_results(results, args.output, args.het_mode)
    
    print("\n" + "=" * 70)
    print("Complete!♡(*´∀｀*)人(*´∀｀*)♡")
    print("=" * 70)


if __name__ == "__main__":
    main()