# Shared structural helper extracted from the original implementation.


LOCAL_LINK = ["EN", "LT", "LM", "LTC", "LMC","LBC"]

LOCAL_CYCLE_F = ["CFML","CFM", "CF", "CFMB", "CFB", "CFL"]

LOCAL_CYCLE_B = ["CFB","CB", "CBL", "CFMB", "CSB"]

LOCAL_CYCLE = ["C","CB","CL","CS","CF","CBL","CFM",
                "CFB","CFL","CSB","CSL","CFMB","CFML",]

SINGLE_BOND = ["S", "eS", "fs", "bS", "s"]

DOUBLE_BOND = ["D", "eD", "fd", "bD", "d"]

TRIPLE_BOND = ["T", "ft", "t"]

CONJU_BOND = ["J", "eJ", "fj", "bJ", "j"]

MULTI_BOND = DOUBLE_BOND + TRIPLE_BOND

BREAK_BOND = ["S", "eS", "bS", "br", "J"]

RING_BOND = ["s","fs","fd","d","r","fr","j","fj","t","ft"]

CHAIN_BOND = ["S","eS","bS","D","eD","bD","T","J","eJ","bJ"]

score_type = [
    [1, 4, 7, 8, 10],
    [2, 3, 5, 6, 9, 11, 12, 13],
    [30, 31, 32],
    [50, 51, 52, 53, 54, 55],
    [60, 61, 62, 63, 64, 65, 66, 67],
]

SINGLE_CONNECT_ATOM_OR_GROUP = ["H", "F", "Cl", "Br"]
