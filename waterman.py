import sys
import numpy as np

def dna_reverse(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    return dna.translate(complement)[::-1]

def dna2num(c):
    mapping = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
    return mapping.get(c, 0)

def dna2num_array(dna):
    char_to_num = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
    return np.array([char_to_num.get(c, 0) for c in dna], dtype=np.int64)

def smith_waterman_multiple(query, ref, match=3, mismatch=-2, gap=-2, min_score=20):
    """使用Smith-Waterman算法查找多个局部比对"""
    m = len(query)
    n = len(ref)
    print(m)
    print(n)
    H = [[0] * (n + 1) for _ in range(m + 1)] # 得分矩阵
    T = [[0] * (n + 1) for _ in range(m + 1)] # 方向矢量
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diag = H[i-1][j-1] + (match if query[i-1] == ref[j-1] else mismatch)
            score_up = H[i-1][j] + gap
            score_left = H[i][j-1] + gap
            H[i][j] = max(0, score_diag, score_up, score_left)
            if H[i][j] == 0:
                T[i][j] = 0
            elif H[i][j] == score_diag:
                T[i][j] = 1
            elif H[i][j] == score_up:
                T[i][j] = 2
            else:
                T[i][j] = 3
    alignments = []
    # print(T)
    # for i in range(1, m + 1):
    #     print(H[i][1:n + 1])
    while True:
        max_score = 0
        max_i, max_j = -1, -1
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if H[i][j] > max_score:
                    max_score = H[i][j]
                    max_i, max_j = i, j
        if max_score < min_score:
            break
        print(max_score, max_i, max_j)
        alignment = traceback(H, T, (max_i, max_j))
        if alignment is not None:
            alignments.append(alignment)
            i, j = max_i, max_j
            while H[i][j] > 0:
                H[i][j] = 0
                if T[i][j] == 1:
                    i -= 1
                    j -= 1
                elif T[i][j] == 2:
                    i -= 1
                else:
                    j -= 1
            break
    return alignments

def traceback(H, T, start_pos):
    """回溯以确定比对区域的起始和结束位置"""
    i, j = start_pos
    query_positions = set()
    ref_positions = set()
    while H[i][j] > 0:
        if T[i][j] == 1:  # 对角线（匹配/失配）
            query_positions.add(i - 1)
            ref_positions.add(j - 1)
            i -= 1
            j -= 1
        elif T[i][j] == 2:  # 向上（查询序列间隙）
            query_positions.add(i - 1)
            i -= 1
        else:  # 向左（目标序列间隙）
            ref_positions.add(j - 1)
            j -= 1
    if not query_positions:
        return None
    query_start = min(query_positions)
    query_end = max(query_positions)
    ref_start = min(ref_positions)
    ref_end = max(ref_positions)
    return (query_start, query_end, ref_start, ref_end)

def find_matching_regions(query, reference, min_score=20):
    """查找查询序列与参考序列的匹配区域，包括直接匹配和倒位匹配"""
    direct_alignments = smith_waterman_multiple(query, reference, min_score=min_score)
    ref_rc = dna_reverse(reference)
    inversion_alignments = smith_waterman_multiple(query, ref_rc, min_score=min_score)
    mapped_inversion_alignments = []
    for align in inversion_alignments:
        qs, qe, rs_rc, re_rc = align
        ref_start = len(reference) - 1 - re_rc
        ref_end = len(reference) - 1 - rs_rc
        mapped_inversion_alignments.append((qs, qe, ref_start, ref_end))
    all_alignments = direct_alignments + mapped_inversion_alignments
    all_alignments.sort(key=lambda x: x[0])
    return all_alignments

# 示例使用
if __name__ == "__main__":
    ref = sys.stdin.readline().strip()
    query = sys.stdin.readline().strip()
    alignments = find_matching_regions(query, ref, min_score=4)  # 低阈值用于测试
    # for align in alignments:
    #     print(align)