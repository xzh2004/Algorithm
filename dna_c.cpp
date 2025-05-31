#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <thread>
#include <future>

const long long MOD = 1000000000007LL;

struct RefSeq {
    int start;
    int end;
    bool reverse;
    
    RefSeq() : start(0), end(0), reverse(false) {}
    RefSeq(int s, int e, bool r) : start(s), end(e), reverse(r) {}
};

struct Trace {
    RefSeq ref_seq;
    int next;
    
    Trace() : next(0) {}
    Trace(const RefSeq& rs, int n) : ref_seq(rs), next(n) {}
};

// DNA reverse complement
std::string dna_reverse(const std::string& dna) {
    std::string result;
    result.reserve(dna.length());
    
    for (int i = dna.length() - 1; i >= 0; i--) {
        switch (dna[i]) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default: result += dna[i]; break;
        }
    }
    return result;
}

// Convert DNA character to number
inline int dna2num(char c) {
    switch (c) {
        case 'A': return 1;
        case 'T': return 2;
        case 'C': return 3;
        case 'G': return 4;
        default: return 0;
    }
}

// Convert DNA string to number array
std::vector<int> dna2num_array(const std::string& dna) {
    std::vector<int> result;
    result.reserve(dna.length());
    for (char c : dna) {
        result.push_back(dna2num(c));
    }
    return result;
}

// Process a chunk of DNA sequence for hashing
std::unordered_map<long long, RefSeq> process_chunk(
    int start, int end, const std::string& dna, bool reverse, int dna_len) {
    
    std::unordered_map<long long, RefSeq> local_map;
    std::string processed_dna = reverse ? dna_reverse(dna) : dna;
    
    for (int i = start; i < end; i++) {
        long long hs = 0;
        if (i % 1000 == 0) printf("Processing chunk: %d to %d\n", i, end);
        for (int j = i + 29; j < dna_len; j++) {
            hs = (hs * 5 + dna2num(processed_dna[j])) % MOD;
            if (local_map.find(hs) == local_map.end()) {
                if (reverse) {
                    local_map[hs] = RefSeq(dna_len - j - 1, dna_len - i - 1, reverse);
                } else {
                    local_map[hs] = RefSeq(i, j, reverse);
                }
            }
        }
    }
    return local_map;
}

// Build hash map for DNA sequence
void hsmake(const std::string& dna, std::unordered_map<long long, RefSeq>& map_, bool reverse) {
    int dna_len = dna.length();
    int num_threads = std::min(8, dna_len / 1000 + 1);
    int chunk_size = dna_len / num_threads;
    
    std::vector<std::future<std::unordered_map<long long, RefSeq>>> futures;
    
    for (int i = 0; i < dna_len; i += chunk_size) {
        int end = std::min(i + chunk_size, dna_len);
        futures.push_back(std::async(std::launch::async, process_chunk, 
                                   i, end, dna, reverse, dna_len));
    }
    
    std::cerr << "Building hash map..." << std::endl;
    for (auto& future : futures) {
        auto local_map = future.get();
        for (const auto& pair : local_map) {
            map_[pair.first] = pair.second;
        }
    }
}

// Search for optimal alignment
std::vector<Trace> search(const std::string& query, 
                         const std::unordered_map<long long, RefSeq>& ref_map) {
    int query_len = query.length();
    std::vector<double> dp(query_len + 1, std::numeric_limits<double>::infinity());
    std::vector<Trace> chain(query_len + 1);
    
    dp[query_len] = 0;
    
    std::vector<int> query_nums = dna2num_array(query);
    
    std::cerr << "Searching sequences..." << std::endl;
    for (int start = query_len - 1; start >= 0; start--) {
        long long hs = 0;
        for (int end = start + 29; end < query_len; end++) {
            hs = (hs * 5 + query_nums[end]) % MOD;
            auto it = ref_map.find(hs);
            if (it != ref_map.end()) {
                const RefSeq& ref_seq = it->second;
                if ((dp[start] > dp[end + 1] + 1) || 
                    (dp[start] == dp[end + 1] + 1 && !ref_seq.reverse)) {
                    dp[start] = dp[end + 1] + 1;
                    chain[start] = Trace(ref_seq, end + 1);
                }
            }
        }
    }
    return chain;
}

// Reconstruct the alignment result
std::vector<std::tuple<int, int, int, int>> reconstruct(
    const std::vector<Trace>& chain, int query_len) {
    
    std::vector<std::tuple<int, int, int, int>> result;
    int pos = 0;
    
    while (pos < query_len) {
        if (chain[pos].next == 0 && pos != query_len) {
            throw std::runtime_error("No answer");
        }
        
        int a = pos;
        int b = chain[pos].next;
        int c = chain[pos].ref_seq.start;
        int d = chain[pos].ref_seq.end;
        
        result.push_back(std::make_tuple(a, b - 1, c, d));
        pos = b;
    }
    
    return result;
}

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    freopen("input3.txt", "r", stdin);
    freopen("output3.txt", "w", stdout);

    std::string ref, query;
    std::getline(std::cin, ref);
    std::getline(std::cin, query);
    
    // Remove any trailing whitespace
    ref.erase(ref.find_last_not_of(" \n\r\t") + 1);
    query.erase(query.find_last_not_of(" \n\r\t") + 1);
    
    std::cerr << "Reference length: " << ref.length() 
              << ", Query length: " << query.length() << std::endl;
    
    std::unordered_map<long long, RefSeq> ref_map;
    std::cerr << "Processing reference sequence..." << std::endl;
    hsmake(ref, ref_map, false);
    hsmake(ref, ref_map, true);
    
    std::cerr << "Processing query sequence..." << std::endl;
    std::vector<Trace> chain = search(query, ref_map);
    auto result = reconstruct(chain, query.length());
    
    // Output results
    std::cout << "[";
    for (size_t i = 0; i < result.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << "(" << std::get<0>(result[i]) << ", " 
                  << std::get<1>(result[i]) << ", "
                  << std::get<2>(result[i]) << ", " 
                  << std::get<3>(result[i]) << ")";
    }
    std::cout << "]" << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cerr << "Total execution time: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}