import sys
from collections import defaultdict

def dna_reverse(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(c, c) for c in reversed(dna))

def dna2num(dna):
    mapping = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
    return mapping.get(dna, 0)

class RefSeq:
    def __init__(self, start, end, reverse):
        self.start = start
        self.end = end
        self.reverse = reverse

MOD = 1_000_000_000_007

def hsmake(dna, map_, reverse):
    dna_len = len(dna)
    if reverse:
        dna = dna_reverse(dna)
    for start in range(dna_len):
        hs = 0
        for end in range(start, dna_len):
            hs = (hs * 5 + dna2num(dna[end])) % MOD
            if hs not in map_:
                if reverse:
                    map_[hs] = RefSeq(dna_len - end - 1, dna_len - start - 1, reverse)
                else:
                    map_[hs] = RefSeq(start, end, reverse)

class Trace:
    def __init__(self, ref_seq, next_):
        self.ref_seq = ref_seq
        self.next = next_

def search(query, ref_map):
    query_len = len(query)
    dp = [float('inf')] * (query_len + 1)
    dp[query_len] = 0
    chain = [None] * (query_len + 1)
    for start in range(query_len - 1, -1, -1):
        hs = 0
        for end in range(start, query_len):
            hs = (hs * 5 + dna2num(query[end])) % MOD
            if hs in ref_map:
                ref_seq = ref_map[hs]
                if (dp[start] > dp[end + 1] + 1) or (dp[start] == dp[end + 1] + 1 and not ref_seq.reverse):
                    dp[start] = dp[end + 1] + 1
                    chain[start] = Trace(ref_seq, end + 1)
    return chain

def reconstruct(chain, query_len):
    result = []
    pos = 0
    while pos < query_len:
        if chain[pos] is None:
            raise ValueError("No answer")
        result.append(chain[pos].ref_seq)
        pos = chain[pos].next
    return result

def main():
    lines = sys.stdin.readlines()
    ref= lines[0].strip()
    query = lines[1].strip()

    ref_map = defaultdict(RefSeq)
    hsmake(ref, ref_map, False)
    hsmake(ref, ref_map, True)

    chain = search(query, ref_map)
    result = reconstruct(chain, len(query))

    for x in result:
        print("POS in REF:", x.start, "Repeat size:", x.end - x.start + 1, "Inverse:", x.reverse)

if __name__ == "__main__":
    main()