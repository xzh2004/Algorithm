import sys
from collections import defaultdict
from tqdm import tqdm
import time
from concurrent.futures import ThreadPoolExecutor
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

class RefSeq:
    __slots__ = ('start', 'end', 'reverse')
    def __init__(self, start, end, reverse):
        self.start = start
        self.end = end
        self.reverse = reverse

MOD = 1_000_000_000_007

def process_chunk(args):
    start, end, dna, reverse, dna_len = args
    local_map = {}
    if reverse:
        dna = dna_reverse(dna)
    for i in range(start, end):
        hs = 0
        for j in range(i, dna_len):
            hs = (hs * 5 + dna2num(dna[j])) % MOD
            if hs not in local_map:
                if reverse:
                    local_map[hs] = RefSeq(dna_len - j - 1, dna_len - i - 1, reverse)
                else:
                    local_map[hs] = RefSeq(i, j, reverse)
    return local_map

def hsmake(dna, map_, reverse):
    dna_len = len(dna)
    num_threads = min(8, dna_len // 1000 + 1)
    chunk_size = dna_len // num_threads
    
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        chunks = [(i, min(i + chunk_size, dna_len), dna, reverse, dna_len) 
                 for i in range(0, dna_len, chunk_size)]
        
        futures = list(tqdm(
            executor.map(process_chunk, chunks),
            total=len(chunks),
            desc="Building hash map",
            file=sys.stderr
        ))
        
        for local_map in futures:
            map_.update(local_map)

class Trace:
    __slots__ = ('ref_seq', 'next')
    def __init__(self, ref_seq, next_):
        self.ref_seq = ref_seq
        self.next = next_

def search(query, ref_map):
    query_len = len(query)
    dp = np.full(query_len + 1, float('inf'), dtype=np.float64)
    dp[query_len] = 0
    chain = [None] * (query_len + 1)
    
    query_nums = dna2num_array(query)
    
    for start in tqdm(range(query_len - 1, -1, -1), desc="Searching sequences", file=sys.stderr):
        hs = 0
        for end in range(start, query_len):
            hs = (hs * 5 + query_nums[end]) % MOD
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
        a, b, c, d = pos, chain[pos].next, chain[pos].ref_seq.start, chain[pos].ref_seq.end
        # while b < query_len and chain[b].next - b < 10 and chain[a].ref_seq.reverse == chain[b].ref_seq.reverse:
        #     d = chain[b].ref_seq.end
        #     b = chain[b].next
        result.append((a, b - 1, c, d))
        # print(f"pos: {pos} -> next: {chain[pos].next} len:{chain[pos].ref_seq.end - chain[pos].ref_seq.start + 1}")
        # if (chain[pos].ref_seq.end - chain[pos].ref_seq.start !=  chain[pos].next - pos):
        #     print(f"pos: {pos} -> next: {chain[pos].next} len:{chain[pos].ref_seq.end - chain[pos].ref_seq.start + 1}")
        # if (b >= query_len): break
        pos = b
    # for x in result:
        
    return result

def main():
    start_time = time.time()
    
    ref = sys.stdin.readline().strip()
    query = sys.stdin.readline().strip()

    print(f"Reference length: {len(ref)}, Query length: {len(query)}", file=sys.stderr)

    ref_map = {}
    print("Processing reference sequence...", file=sys.stderr)
    hsmake(ref, ref_map, False)
    hsmake(ref, ref_map, True)

    print("Processing query sequence...", file=sys.stderr)
    chain = search(query, ref_map)
    result = reconstruct(chain, len(query))

    # pos = 0
    # output = []
    # for x in result:
    #     output.append((pos, pos + x.end - x.start, x.start, x.end))
    #     pos += x.end - x.start + 1
    # print(f"pos: {pos}")
    end_time = time.time()
    print(result)
    print(f"Total execution time: {end_time - start_time:.2f} seconds", file=sys.stderr)

if __name__ == "__main__":
    main()