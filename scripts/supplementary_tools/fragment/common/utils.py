# Shared structural helper extracted from the original implementation.
from typing import List

def combine_arr(a: List[list]) -> List[list]:
    """Input a is a 2-D array. Merge all arrays which shares at least 1 common element
    e.g, [[1,2], [2,3], [4,5], [6,7], [8,9], [7,8]]
    returns [[1,2,3], [4,5], [6,7,8,9]]
    """
    b = len(a)
    for i in range(b):
        if a[i] == []:
            continue
        for j in range(b):
            if i == j:
                break
            if a[j] == []:
                continue
            merged = list(set(a[i] + a[j]))
            total_len = len(a[j]) + len(a[i])
            if len(merged) < total_len:  # there is intersection
                a[i], a[j] = merged, []
    return [ele for ele in a if ele != []]
