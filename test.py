import sys

def compare_substrings_by_hash(ref, query, query_start, query_end, ref_start, ref_end):
    """比较两个子串是否相同，通过比较哈希值"""
    # 提取子串
    query_sub = query[query_start:query_end+1]
    ref_sub = ref[ref_start:ref_end+1]
    
    # 计算哈希值
    query_hash = hash(query_sub)
    ref_hash = hash(ref[ref_start:ref_end+1])
    
    if query_hash != ref_hash:
        print(f"len:{len(query_sub)}")
        print(query_sub, ref_sub)

    # 比较哈希值
    if query_hash == ref_hash:
        # 哈希值相同，再比较实际内容以确保正确性
        return query_sub == ref_sub
    return False

def main():
    # 读取输入
    ref = sys.stdin.readline().strip()
    query = sys.stdin.readline().strip()
    print(query[-50:])
    print(ref[-50:])
    # 处理每一组比较
    # for line in sys.stdin:
    #     line = line.strip()
    #     if not line:
    #         continue
            
    #     # 解析坐标
    #     coords = line.strip('()').split(',')
    #     query_start = int(coords[0].strip())
    #     query_end = int(coords[1].strip())
    #     ref_start = int(coords[2].strip())
    #     ref_end = int(coords[3].strip())
        
    #     # 比较并输出结果
    #     is_same = compare_substrings_by_hash(ref, query, query_start, query_end, ref_start, ref_end)
    #     if (is_same == 0): print(f"{query_start}-{query_end} vs {ref_start}-{ref_end}: {'yes' if is_same else 'no'}")

if __name__ == "__main__":
    main()
