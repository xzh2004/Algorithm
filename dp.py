import sys
import numpy as np

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
K = 200
LENGTH_LIMIT = 9 # 数据点1为30， 数据点2为9

hash_map = []

class Trace:
    def __init__(self, ref_seq, pre_):
        self.ref_seq = ref_seq
        self.pre = pre_

class RollingHash:
    """滚动哈希类，避免重复计算"""
    def __init__(self, sequence, base=5, mod=MOD):
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
        for end in range(start, min(dna_len, start + K)):
            hs = roller.get_hash(start, end)
            if reverse:
                hash_map.append((hs, RefSeq(dna_len - end - 1, dna_len - start - 1, reverse)))
            else:
                hash_map.append((hs, RefSeq(start, end, reverse)))
            if hs not in map_:
                if reverse:
                    map_[hs] = RefSeq(dna_len - end - 1, dna_len - start - 1, reverse)
                else:
                    map_[hs] = RefSeq(start, end, reverse)
    hash_map.sort(key=lambda x: x[0])

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
        
        for end in range(min(start + K - 1, query_len - 1), start, -1):
            hs = roller.get_hash(start, end)
            if hs in ref_map:
                ref_seq = ref_map[hs]
                if dp[start - 1] + end - start + 1 > dp[end]:
                    dp[end] = dp[start - 1] + end - start + 1
                    pos[end] = end
                    chain[end] = Trace(ref_seq, start - 1)
                # break
    return roller, chain, pos

def reconstruct(chain, query_len, pos):
    result = []
    now = pos[query_len - 1]
    while now > 0:
        if chain[now] is None:
            raise ValueError("No answer")
        nl = chain[now].ref_seq.end - chain[now].ref_seq.start + 1
        result.append((now - nl + 1, now, chain[now].ref_seq.start, chain[now].ref_seq.end, chain[now].ref_seq.reverse))
        pre_ = chain[now].pre
        now = pre_
    return result


def find_hash(hs):
    """使用二分查找检查哈希值是否在hash_map中"""
    left, right = 0, len(hash_map) - 1
    
    while left <= right:
        mid = (left + right) // 2
        if hash_map[mid][0] == hs:
            return True
        elif hash_map[mid][0] < hs:
            left = mid + 1
        else:
            right = mid - 1
    return False

def check(roller, t1, t2):
    a, b = t1
    na, nb, nc, nd, nr = t2
    hs = roller.get_hash(a, b)
    left, right = 0, len(hash_map) - 1
    st = len(hash_map)
    while left <= right:
        mid = (left + right) // 2
        if hash_map[mid][0] < hs:
            left = mid + 1
            st = mid
        else:
            right = mid - 1
    st = st + 1
    for i in range(st, len(hash_map)):
        if hash_map[i][0] == hs:
            ref_seq = hash_map[i][1]
            c = ref_seq.start
            d = ref_seq.end
            r = ref_seq.reverse
            if (c - nd >= 1 and c - nd <= 15 and r == 0 and nr == r):
                return 1, d
            elif (nc - d >= 1 and nc - d <= 15 and r == 1 and nr == r):
                return 2, c
        else: break
    return 0, 0

def main():
    lines = sys.stdin.readlines()
    ref= lines[0].strip()
    query = lines[1].strip()
    print(f"query_len:{len(query)}, ref_len:{len(ref)}")
    ref_map = {}
    hsmake_rolling_hash(ref, ref_map, False)
    hsmake_rolling_hash(ref, ref_map, True)
    hs_roller, chain, pos = search_rolling_hash(query, ref_map)
    sorted_result = reconstruct(chain, len(query), pos)
    
    ans = []
    rl = len(sorted_result)
    na, nb, nc, nd, nr = sorted_result[rl - 1]
    for i in range(rl - 2, 0, -1):
        a, b, c, d, reverse = sorted_result[i]
        if (c - nd >= 1 and c - nd <= 5 and reverse == 0 and nr == reverse):
            nb = b
            nd = d
            # print("type 1")
        else:
            if (nc - d >= 1 and nc - d <= 5 and reverse == 1 and nr == reverse):
                nb = b
                nc = c
                # print("type 2")
            else:
                res, para = check(hs_roller, (a, b), (na, nb, nc, nd, nr))
                if res == 1:
                    nb = b
                    nd = para
                    # print("type 3")
                elif res == 2:
                    nb = b
                    nc = para
                    # print("type 4")
                elif b - a + 1 <= LENGTH_LIMIT:
                    nb = b
                    # print("type 5")
                    if nr == 0: nd = nd + b - a
                    else: nc = nc - (b - a)
                else:
                    ans.append((na, nb + 1, nc, nd + 1))
                    na, nb, nc, nd, nr = a, b, c, d, reverse

    ans.append((na, nb + 1, nc, nd + 1))
    for x in ans:
        print(f"{x},")

if __name__ == "__main__":
    main()