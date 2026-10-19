from io import StringIO

import numpy as np
import pandas as pd

alphabet = "ACGT"
transitions = ["A>G", "C>T", "G>A", "T>C"]
transversions = ["A>C", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>G"]
color_mapping12 = {
    "C>A": "deepskyblue",
    "G>T": "deepskyblue",
    "C>G": "black",
    "G>C": "black",
    "C>T": "red",
    "G>A": "red",
    "T>A": "silver",
    "A>T": "silver",
    "T>C": "yellowgreen",
    "A>G": "yellowgreen",
    "T>G": "pink",
    "A>C": "pink",
}
amino_acid_codes = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}

grantham_matrix_str = """
FIRST	R	L	P	T	A	V	G	I	F	Y	C	H	Q	N	K	D	E	M	W
S	110	145	74	58	99	124	56	142	155	144	112	89	68	46	121	65	80	135	177
R	0	102	103	71	112	96	125	97	97	77	180	29	43	86	26	96	54	91	101
L	0	0	98	92	96	32	138	5	22	36	198	99	113	153	107	172	138	15	61
P	0	0	0	38	27	68	42	95	114	110	169	77	76	91	103	108	93	87	147
T	0	0	0	0	58	69	59	89	103	92	149	47	42	65	78	85	65	81	128
A	0	0	0	0	0	64	60	94	113	112	195	86	91	111	106	126	107	84	148
V	0	0	0	0	0	0	109	29	50	55	192	84	96	133	97	152	121	21	88
G	0	0	0	0	0	0	0	135	153	147	159	98	87	80	127	94	98	127	184
I	0	0	0	0	0	0	0	0	21	33	198	94	109	149	102	168	134	10	61
F	0	0	0	0	0	0	0	0	0	22	205	100	116	158	102	177	140	28	40
Y	0	0	0	0	0	0	0	0	0	0	194	83	99	143	85	160	122	36	37
C	0	0	0	0	0	0	0	0	0	0	0	174	154	139	202	154	170	196	215
H	0	0	0	0	0	0	0	0	0	0	0	0	24	68	32	81	40	87	115
Q	0	0	0	0	0	0	0	0	0	0	0	0	0	46	53	61	29	101	130
N	0	0	0	0	0	0	0	0	0	0	0	0	0	0	94	23	42	142	174
K	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	101	56	95	110
D	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	160	181
E	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	126	152
M	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	67
"""

grantham = pd.read_csv(StringIO(grantham_matrix_str), sep="\t", index_col=0).replace(0, np.nan)
grantham.index.name = "aa1"
grantham_long = grantham.melt(
    ignore_index=False, var_name="aa2", value_name="grantham_distance"
).dropna().reset_index()
grantham_long["aa1"] = grantham_long["aa1"].map(amino_acid_codes)
grantham_long["aa2"] = grantham_long["aa2"].map(amino_acid_codes)
grantham_long = pd.concat(
    [grantham_long, grantham_long.rename(columns={"aa1": "aa2", "aa2": "aa1"})],
    ignore_index=True,
)
del grantham
