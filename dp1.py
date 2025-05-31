# data1
import sys
import numpy as np
from collections import defaultdict

def dna_reverse(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    return dna.translate(complement)[::-1]

def dna2num(dna):
    mapping = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
    return mapping.get(dna, 0)

def dna_reverse(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    return dna.translate(complement)[::-1]

def dna2num_vectorized(dna_str):
    """向量化的DNA到数字转换"""
    mapping = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
    return np.array([mapping.get(c, 0) for c in dna_str], dtype=np.int8)

class RefSeq:
    def __init__(self, start, end, reverse):
        self.start = start
        self.end = end
        self.reverse = reverse

MOD = 1_000_000_000_007
INV = 4_000_000_000_03

def hsmake(dna, map_, reverse):
    dna_len = len(dna)
    if reverse:
        dna = dna_reverse(dna)
    for start in range(dna_len):
        hs = 0
        # for i in range(start, min(dna_len, start + 28)):
        #     hs = (hs * 5 + dna2num(dna[i])) % MOD
        for end in range(start, min(dna_len, start)):
            hs = (hs * 5 + dna2num(dna[end])) % MOD
            if (start == 1 and end == 60):
                print(f"start:{start} end:{end} hs:{hs}")
            if hs not in map_:
                if reverse:
                    map_[hs] = RefSeq(dna_len - end - 1, dna_len - start - 1, reverse)
                else:
                    map_[hs] = RefSeq(start, end, reverse)

class Trace:
    def __init__(self, ref_seq, pre_):
        self.ref_seq = ref_seq
        self.pre = pre_

def search(query, ref_map):
    query_len = len(query)
    dp = [0] * (query_len + 5)
    pos  = [0] * (query_len + 5)
    pow5 = [0] * (query_len + 5)
    pow5[0] = 1
    chain = [None] * (query_len + 5)
    for i in range(1, query_len + 5):
        pow5[i] = (pow5[i - 1] * 5) % MOD
    for start in range(query_len):
        hs = 0
        if start > 0 and dp[start] < dp[start - 1]:
            pos[start] = pos[start - 1]
            dp[start] = dp[start - 1]
        # print(f"start:{start} pos:{pos[start]}")
        for i in range(start, query_len):
            hs = (hs * 5 + dna2num(query[i])) % MOD
        for end in range(query_len - 1, start, -1):
            # hs = (hs * 5 + dna2num(query[end])) % MOD
            # print(f"query: {start} - {end}")
            # if (start == 1 and end == 60):
            #     print(f"start:{start} end:{end} hs:{hs}")
            if hs in ref_map:
                ref_seq = ref_map[hs]
                # print(f"query_matched: {start} - {end}")
                # if (dp[start] > dp[end + 1] + 1) or (dp[start] == dp[end + 1] + 1 and not ref_seq.reverse):
                #     dp[start] = dp[end + 1] + 1
                #     chain[start] = Trace(ref_seq, end + 1)
                if (dp[start - 1] + end - start + 1 > dp[end]):
                    # print(f"query_dp: {start} - {end}")
                    dp[end] = dp[start - 1] + end - start + 1
                    pos[end] = end
                    chain[end] = Trace(ref_seq, start - 1)
                break
            hs = (hs - dna2num(query[end]) + MOD) % MOD * INV % MOD
    return chain, pos

class RollingHash:
    """滚动哈希类，避免重复计算"""
    def __init__(self, sequence, base=5, mod=1_000_000_000_007):
        self.sequence = sequence
        self.base = base
        self.mod = mod
        self.mapping = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
        
        # 预计算base的幂次
        self.powers = [1]
        for i in range(1, len(sequence) + 1):
            self.powers.append((self.powers[-1] * base) % mod)
        
        # 预计算前缀哈希
        self.prefix_hash = [0]
        for char in sequence:
            val = self.mapping.get(char, 0)
            new_hash = (self.prefix_hash[-1] * base + val) % mod
            self.prefix_hash.append(new_hash)
    
    def get_hash(self, start, end):
        """获取子串[start:end+1]的哈希值"""
        if start == 0:
            return self.prefix_hash[end + 1]
        
        # hash[start:end+1] = prefix[end+1] - prefix[start] * base^(end-start+1)
        length = end - start + 1
        result = (self.prefix_hash[end + 1] - 
                 self.prefix_hash[start] * self.powers[length]) % self.mod
        return (result + self.mod) % self.mod  # 确保非负

def hsmake_rolling_hash(dna, map_, reverse):
    """使用滚动哈希的版本"""
    if reverse:
        dna = dna_reverse(dna)
    
    roller = RollingHash(dna)
    dna_len = len(dna)
    
    for start in range(dna_len):
        for end in range(start, min(dna_len, start + 200)):
            hs = roller.get_hash(start, end)
            if (start == 1 and end == 60):
                print(f"start:{start} end:{end} hs:{hs}")
            if hs not in map_:
                if reverse:
                    map_[hs] = RefSeq(dna_len - end - 1, dna_len - start - 1, reverse)
                else:
                    map_[hs] = RefSeq(start, end, reverse)

def search_rolling_hash(query, ref_map):
    """使用滚动哈希的搜索"""
    query_len = len(query)
    roller = RollingHash(query)
    
    dp = [0] * (query_len + 5)
    pos = [0] * (query_len + 5)
    chain = [None] * (query_len + 5)
    
    for start in range(query_len):
        if start > 0 and dp[start] < dp[start - 1]:
            pos[start] = pos[start - 1]
            dp[start] = dp[start - 1]
        
        for end in range(min(start + 199, query_len - 1), start, -1):
            hs = roller.get_hash(start, end)
            if (start == 1 and end == 60):
                print(f"start:{start} end:{end} hs:{hs}")
            if hs in ref_map:
                ref_seq = ref_map[hs]
                # print(f"query_matched: {start} - {end}")
                if dp[start - 1] + end - start + 1 > dp[end]:
                    dp[end] = dp[start - 1] + end - start + 1
                    pos[end] = end
                    chain[end] = Trace(ref_seq, start - 1)
                # break
    
    return chain, pos

def reconstruct(chain, query_len, pos):
    result = []
    now = pos[query_len - 1]
    # print(f"now:{now}")
    while now > 0:
        # print(now)
        if chain[now] is None:
            raise ValueError("No answer")
        nl = chain[now].ref_seq.end - chain[now].ref_seq.start + 1
        result.append((now - nl + 1, now, chain[now].ref_seq.start, chain[now].ref_seq.end, chain[now].ref_seq.reverse))
        pre_ = chain[now].pre
        print(f"len: {nl}", now - nl + 1, now, chain[now].ref_seq.start, chain[now].ref_seq.end)
        now = pre_
    return result
        

def main():
    lines = sys.stdin.readlines()
    ref= lines[0].strip()
    query = lines[1].strip()
    print(f"query_len:{len(query)}, ref_len:{len(ref)}")
    ref_map = {}
    hsmake_rolling_hash(ref, ref_map, False)
    hsmake_rolling_hash(ref, ref_map, True)
    chain, pos = search_rolling_hash(query, ref_map)
    result = reconstruct(chain, len(query), pos)
    # hsmake(ref, ref_map, False)
    # hsmake(ref, ref_map, True)

    # chain, pos = search(query, ref_map)
    # result = reconstruct(chain, len(query), pos)

    ans = []
    rl = len(result)
    na, nb, nc, nd, nr = result[rl - 1]
    for i in range(rl - 2, 0, -1):
        a, b, c, d, reverse = result[i]
        if (c - nd >= 1 and c - nd <= 5 and reverse == 0 and nr == reverse):
            nb = b
            nd = d
        else:
            if (nc - d >= 1 and nc - d <= 5 and reverse == 1 and nr == reverse):
                nb = b
                nc = c
            else:
                if b - a + 1 <= 9:
                    nb = b
                    if nr == 0: nd = nd + b - a
                    else: nc = nc - (b - a)
                else:
                    ans.append((na, nb + 1, nc, nd + 1))
                    na, nb, nc, nd, nr = a, b, c, d, reverse
        # print(f"{a}, {b}, {c}, {d} -> {na}, {nb}, {nc}, {nd}")
    ans.append((na, nb + 1, nc, nd + 1))
    for x in ans:
        print(f"{x},")

if __name__ == "__main__":
    main()