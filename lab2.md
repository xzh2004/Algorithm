# 算法设计与分析 lab2 实验文档

姓名：项正豪  学号：23307130156

github仓库地址：https://github.com/xzh2004/Algorithm

## 算法思路

先对ref计算hash值，然后在query中寻找能与ref匹配的完全相同的hash值。以部分能够完美匹配的query的子串作为锚点（记录为chain），向左右延伸。因为不是所有完美匹配的子串都能作为锚点，所以采用动态规划的策略确认锚点。设状态dp[i]表示匹配到query中的i位置为止，匹配上的最大长度，pos[i]记录转移到这个状态的上一个状态的区间右端点。

状态转移方程为：dp[end] = max(dp[end - 1],  dp[start - 1] + end - start + 1)，其中start和end为每一个能够完美匹配的query的子串。同时也会对pos和chain作出相应的更新。这一状态转移方程借鉴了**最短路中路径松弛的思想**，所以可以认为是使用了图相关的算法

然后通过pos和chain进行重构，得到一个全部是完美匹配的，且能实现对query串最长覆盖的区间序列，类似于下图。

![dna_alignment](./dna_alignment.png)

图中橙色的为此时得到的区间序列，可以看到匹配效果还不错，但有很多都是断开的，所以我们还需考虑对这些橙色序列进行合并。

合并序列的逻辑用以下伪代码描述：

```python
# 定义变量：
# (a,b): 当前匹配在查询序列的位置 [a,b]
# (c,d): 当前匹配在参考序列的位置 [c,d]
# (na,nb): 已合并的查询序列位置 [na,nb]
# (nc,nd): 已合并的参考序列位置 [nc,nd]
# reverse: 当前匹配的方向
# nr: 已合并匹配的方向

# 合并条件分类：

# 类型1：正向匹配，参考序列位置间隔在[1,5]之间
if (c - nd ∈ [1,5]) ∧ (reverse = 0) ∧ (nr = reverse):
    更新: nb = b
          nd = d

# 类型2：反向匹配，参考序列位置间隔在[1,5]之间
elif (nc - d ∈ [1,5]) ∧ (reverse = 1) ∧ (nr = reverse):
    更新: nb = b
          nc = c

# 类型3：通过check函数检查，返回类型1

elif check() = (1, para):
    更新: nb = b
          nd = para

# 类型4：通过check函数检查，返回类型2
elif check() = (2, para):
    更新: nb = b
          nc = para

# 类型5：短序列处理（长度 ≤ LENGTH_LIMIT）
elif (b - a + 1) ≤ LENGTH_LIMIT:
    更新: nb = b
    if nr = 0:
        nd = nd + (b - a)
    else:
        nc = nc - (b - a)

# 无法合并：添加当前合并结果，开始新的合并
else:
    添加 (na, nb+1, nc, nd+1) 到结果集
    重置: (na,nb,nc,nd,nr) = (a,b,c,d,reverse)
```

其中check函数的作用为：

 check函数的作用为防止错误的完美匹配发生

> e.g. 
>
> ref(1, 5) AGCTA
>
> ref(101,105) AGCTA
>
> query(11, 15) AGCTA

在这个例子中，query(11,15)可能与ref(1,5)匹配，但不会与ref(101,105)匹配，但实际上我们需要它与ref(101, 105)匹配，因为这样才能通过判断ref的区间端点决定是否合并（如类型1和2中做的那样）。为了解决这个问题，我手写实现了一个可重集hash_map，在这个可重集中一个hash值可以对应多个ref序列，并且支持二分查找。

以上，合并后的序列即为最终的答案序列

## 算法伪代码

```python
# 算法：DNA序列比对与合并算法
# 目标：在参考序列中找到查询序列的最佳匹配位置，支持正向和反向匹配

# ================ 数据结构定义 ================
class RefSeq:
    start: int      # 序列起始位置
    end: int        # 序列结束位置
    reverse: bool   # 是否反向匹配

class Trace:
    ref_seq: RefSeq  # 参考序列信息
    pre: int         # 前一个匹配位置

class RollingHash:
    # 用于高效计算DNA序列子串的哈希值
    sequence: str    # DNA序列
    base: int = 5    # 哈希基数
    mod: int = 1e12+7  # 模数
    
    # 预计算数组
    powers: List[int]      # base的幂次数组
    prefix_hash: List[int] # 前缀哈希数组
    
    function __init__(sequence):
        # 初始化DNA到数字的映射
        mapping = {'A':1, 'T':2, 'C':3, 'G':4}
        
        # 预计算base的幂次
        powers[0] = 1
        for i in range(1, len(sequence)+1):
            powers[i] = (powers[i-1] * base) % mod
        
        # 预计算前缀哈希
        prefix_hash[0] = 0
        for char in sequence:
            val = mapping[char]
            prefix_hash[i+1] = (prefix_hash[i] * base + val) % mod
    
    function get_hash(start, end):
        # 计算子串[start:end+1]的哈希值
        length = end - start + 1
        if start == 0:
            return prefix_hash[end+1]
        
        result = (prefix_hash[end+1] - 
                 prefix_hash[start] * powers[length]) % mod
        return (result + mod) % mod  # 确保非负

# ================ 辅助函数 ================
function dna_reverse(dna):
    # 计算DNA序列的反向互补序列
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(complement[c] for c in dna[::-1])

function dna2num_vectorized(dna_str):
    # 将DNA序列转换为数字向量
    mapping = {'A':1, 'T':2, 'C':3, 'G':4}
    return [mapping[c] for c in dna_str]

# ================ 核心算法实现 ================
function hsmake_rolling_hash(dna, ref_map, reverse):
    # 预处理参考序列，生成哈希表
    if reverse:
        dna = dna_reverse(dna)
    
    roller = RollingHash(dna)
    dna_len = len(dna)
    
    # 处理所有可能的子串
    for start in range(dna_len):
        for end in range(start, min(dna_len, start + 200)):
            hs = roller.get_hash(start, end)
            
            # 存储哈希值和对应的序列信息
            if reverse:
                seq_info = RefSeq(dna_len-end-1, dna_len-start-1, reverse)
            else:
                seq_info = RefSeq(start, end, reverse)
            
            hash_map.append((hs, seq_info))
            if hs not in ref_map:
                ref_map[hs] = seq_info
    
    # 对哈希表排序，便于后续二分查找
    hash_map.sort(key=lambda x: x[0])

function search_rolling_hash(query, ref_map):
    # 使用动态规划进行序列比对
    query_len = len(query)
    roller = RollingHash(query)
    
    # 初始化动态规划数组
    dp = [0] * (query_len + 5)    # 存储最大匹配长度
    pos = [0] * (query_len + 5)   # 存储匹配位置
    chain = [None] * (query_len + 5)  # 存储回溯信息
    
    # 动态规划过程
    for start in range(query_len):
        # 处理重叠情况
        if start > 0 and dp[start] < dp[start-1]:
            pos[start] = pos[start-1]
            dp[start] = dp[start-1]
        
        # 搜索可能的匹配
        for end in range(min(start+199, query_len-1), start, -1):
            hs = roller.get_hash(start, end)
            if hs in ref_map:
                ref_seq = ref_map[hs]
                new_length = dp[start-1] + (end-start+1)
                
                # 更新动态规划状态
                if new_length > dp[end]:
                    dp[end] = new_length
                    pos[end] = end
                    chain[end] = Trace(ref_seq, start-1)
    
    return roller, chain, pos

function reconstruct(chain, query_len, pos):
    # 重构比对结果
    result = []
    current_pos = pos[query_len-1]
    
    while current_pos > 0:
        if chain[current_pos] is None:
            raise ValueError("无法找到有效匹配")
        
        # 获取当前匹配信息
        ref_seq = chain[current_pos].ref_seq
        match_length = ref_seq.end - ref_seq.start + 1
        
        # 添加匹配信息到结果
        result.append((
            current_pos - match_length + 1,  # 查询序列起始位置
            current_pos,                     # 查询序列结束位置
            ref_seq.start,                   # 参考序列起始位置
            ref_seq.end,                     # 参考序列结束位置
            ref_seq.reverse                  # 是否反向匹配
        ))
        
        current_pos = chain[current_pos].pre
    
    return result

function merge_matches(sorted_result):
    # 合并相邻的匹配
    ans = []
    if not sorted_result:
        return ans
    
    # 初始化第一个匹配
    na, nb, nc, nd, nr = sorted_result[-1]
    
    # 处理后续匹配
    for i in range(len(sorted_result)-2, -1, -1):
        a, b, c, d, reverse = sorted_result[i]
        
        # 检查是否可以合并
        if (可以合并条件):
            # 更新当前匹配范围
            nb = b
            if reverse == 0:
                nd = d
            else:
                nc = c
        else:
            # 添加新的匹配
            ans.append((na, nb+1, nc, nd+1))
            na, nb, nc, nd, nr = a, b, c, d, reverse
    
    # 添加最后一个匹配
    ans.append((na, nb+1, nc, nd+1))
    return ans

# ================ 主函数 ================
function main():
    # 读取输入
    ref_seq = 读取参考序列
    query_seq = 读取查询序列
    
    # 初始化数据结构
    ref_map = {}
    hash_map = []
    
    # 预处理参考序列
    hsmake_rolling_hash(ref_seq, ref_map, False)  # 正向处理
    hsmake_rolling_hash(ref_seq, ref_map, True)   # 反向处理
    
    # 执行序列比对
    roller, chain, pos = search_rolling_hash(query_seq, ref_map)
    
    # 重构结果
    sorted_result = reconstruct(chain, len(query_seq), pos)
    
    # 合并相邻匹配
    final_result = merge_matches(sorted_result)
    
    # 输出结果
    for match in final_result:
        输出匹配信息
```

## 时间复杂度

设ref的长度为n，query的长度为m

预处理ref hash（即hsmake_rolling_hash函数）的复杂度为 `O(Kn)`，其中`K`为常数。如果匹配所有区间的hash值，时间复杂度按理说应为`O(nm)`。但考虑到很多过长的区间不可能有完美匹配的情况，所以我们将区间长度限制在一个常数K（本代码中K=200）内，这样可以将时间复杂度压缩为`O(Kn)`，即**线性复杂度**

在query序列中寻找锚点（即search_rolling_hash函数）的复杂度为`O(Km)`, 其中`K`为常数

重构结果（即reconstruct函数）和合并相邻匹配（即merge_matches）的复杂度为`O(m)`

手写实现可重集的复杂度为`O(Kmlog(Km))`，但实际存储的序列个数远没有Km个，所以可近似认为复杂度小于`O(K(n+m))`

综上，总复杂度为`O(K(n+m))`，**线性复杂度**

## 空间复杂度

空间开销主要在存储hash值的字典里，所以空间总复杂度为`O(Kn)`

## 运行结果

### 数据点1

```
[(0, 5172, 0, 5172),
(5172, 5218, 24617, 24663),
(5218, 6678, 5218, 6678),
(6678, 23175, 6653, 23151),
(23175, 29837, 23176, 29837)]
```

得分为29821

### 数据点2

```
[(0, 301, 0, 301),
(301, 400, 401, 500),
(400, 501, 499, 600),
(507, 694, 607, 793),
(694, 788, 694, 788),
(788, 900, 688, 800),
(900, 1000, 700, 800),
(1000, 1200, 700, 900),
(1200, 1308, 894, 1000),
(1308, 1402, 908, 1000),
(1402, 1500, 402, 500),
(1501, 1602, 1001, 1102),
(1602, 1700, 1302, 1400),
(1700, 1802, 1200, 1302),
(1810, 1900, 1110, 1200),
(1900, 1998, 1400, 1498),
(2299, 2500, 1499, 1700)]
```

得分为2092

## 运行说明

程序中LENGTH_LIMIT变量的作用是区间长度小于等于LENGTH_LIMIT的会直接与相邻区间合并

对于第一个数据点，LENGTH_LIMIT的值为30

对于第二个数据点，LENGTH_LIMIT的值为9

运行两个数据点时仅需要修改该值即可，其他的代码都是相同的

## 实验讨论

1. 在考虑区间合并的时候，我想出了一个较为完善的**最长路算法**。将完美匹配的ref和query的左端点和右端点建边，设置边的权值为匹配的长度。相邻的query和ref序列内可以建边（比如端点距离在5以内），这时边的权值为0。

   类似于下图，- 和 | 均表示建的图的边。

   ```
   ref:    ATC  - AAAGTC - T___ - ATGC -
   		       |    |          |  |
   query:  GCT  - AAAGTC - AAAA - ATGC - 
   ```

   

   考虑到目前的代码正确率尚可，以及那个最长路算法代码量比较大，所以没有实现完整的最长路匹配算法，但已经实现的代码中借鉴了这一想法，本实验采用的**动态规划算法类似于最短路中的路径松弛**，也可以说是使用了一种**图的算法**

   

2. 我认为目前代码的短板在于对于区间端点的处理还不够细致，未来可以再完善区间合并的逻辑，以提高得分。

