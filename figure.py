import re
from collections import Counter

def parse_data_file(filename):
    """解析数据文件，提取(query_start, query_end, ref_start, ref_end)"""
    data = []
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
        # 使用正则表达式匹配所有的元组
        pattern = r'\((\d+),\s*(\d+),\s*(\d+),\s*(\d+)\)'
        matches = re.findall(pattern, content)
        
        for match in matches:
            query_start, query_end, ref_start, ref_end = map(int, match)
            data.append((query_start, query_end, ref_start, ref_end))
    
    return data

def create_heatmap_matrix(data, grid_size=50):
    """创建热力图矩阵"""
    query_starts = [item[0] for item in data]
    ref_starts = [item[2] for item in data]
    
    # 获取范围
    min_query, max_query = min(query_starts), max(query_starts)
    min_ref, max_ref = min(ref_starts), max(ref_starts)
    
    # 创建网格
    query_step = (max_query - min_query) // grid_size if grid_size > 0 else 1
    ref_step = (max_ref - min_ref) // grid_size if grid_size > 0 else 1
    
    # 初始化矩阵
    matrix = [[0 for _ in range(grid_size + 1)] for _ in range(grid_size + 1)]
    
    # 填充矩阵
    for query_start, ref_start in zip(query_starts, ref_starts):
        query_idx = min((query_start - min_query) // query_step, grid_size) if query_step > 0 else 0
        ref_idx = min((ref_start - min_ref) // ref_step, grid_size) if ref_step > 0 else 0
        matrix[query_idx][ref_idx] += 1
    
    return matrix, (min_query, max_query, query_step), (min_ref, max_ref, ref_step)

def print_ascii_heatmap(matrix, query_info, ref_info, max_display_size=30):
    """打印ASCII热力图"""
    min_query, max_query, query_step = query_info
    min_ref, max_ref, ref_step = ref_info
    
    print("\n🔥 Query Start vs Reference Start 热力图")
    print("=" * 80)
    print("横坐标: Reference Start位置, 纵坐标: Query Start位置")
    print("数字表示该位置的匹配数量")
    print("-" * 80)
    
    # 获取最大值用于标准化
    max_val = max(max(row) for row in matrix)
    if max_val == 0:
        max_val = 1
    
    # 缩放矩阵以适应显示
    display_size = min(len(matrix), max_display_size)
    step = len(matrix) // display_size if display_size > 0 else 1
    
    # 打印列标题
    print("Query\\Ref", end="")
    for j in range(0, len(matrix[0]), max(1, len(matrix[0]) // 15)):
        ref_pos = min_ref + j * ref_step
        print(f"{ref_pos:>6}", end="")
    print()
    
    # 打印行
    for i in range(0, len(matrix), max(1, len(matrix) // 20)):
        query_pos = min_query + i * query_step
        print(f"{query_pos:>7}", end=" ")
        
        for j in range(0, len(matrix[0]), max(1, len(matrix[0]) // 15)):
            val = matrix[i][j]
            if val == 0:
                print("  .", end="   ")
            else:
                # 使用不同字符表示不同强度
                intensity = val * 9 // max_val
                chars = ['.', '░', '▒', '▓', '█', '🔥']
                char_idx = min(intensity * len(chars) // 10, len(chars) - 1)
                print(f"{val:>3}", end="   ")
        print()

def create_detailed_heatmap(data):
    """创建详细的热力图分析"""
    query_starts = [item[0] for item in data]
    query_ends = [item[1] for item in data]
    ref_starts = [item[2] for item in data]
    ref_ends = [item[3] for item in data]
    
    print("\n📊 详细匹配分析")
    print("=" * 60)
    
    # 按记录索引分析
    print("\n🎯 按记录索引的匹配模式:")
    print("Index | Query Range    | Ref Range      | Match Type")
    print("-" * 55)
    
    for i in range(min(20, len(data))):  # 只显示前20条
        q_start, q_end, r_start, r_end = data[i]
        match_type = ""
        if q_start == r_start and q_end == r_end:
            match_type = "完全匹配"
        elif q_start == r_start:
            match_type = "起始匹配"
        elif q_end == r_end:
            match_type = "结束匹配"
        else:
            match_type = "不匹配"
        
        print(f"{i:5d} | {q_start:4d}-{q_end:4d}     | {r_start:4d}-{r_end:4d}     | {match_type}")
    
    if len(data) > 20:
        print(f"... 还有 {len(data) - 20} 条记录")
    
    # 创建位置对应关系热力图
    print(f"\n🔥 位置对应关系热力图 (Query Start vs Ref Start)")
    print("=" * 70)
    
    # 统计每个(query_start, ref_start)对的出现次数
    position_pairs = {}
    for q_start, r_start in zip(query_starts, ref_starts):
        key = (q_start, r_start)
        position_pairs[key] = position_pairs.get(key, 0) + 1
    
    # 按出现频率排序
    sorted_pairs = sorted(position_pairs.items(), key=lambda x: x[1], reverse=True)
    
    print("Query Start | Ref Start | 出现次数 | 可视化")
    print("-" * 50)
    
    max_count = max(position_pairs.values()) if position_pairs else 1
    for (q_start, r_start), count in sorted_pairs[:25]:  # 显示前25个最频繁的
        bar_length = count * 20 // max_count
        bar = "█" * bar_length
        print(f"{q_start:>10} | {r_start:>8} | {count:>6} | {bar}")

def analyze_matching_patterns(data):
    """分析匹配模式"""
    print(f"\n🔍 匹配模式深度分析")
    print("=" * 60)
    
    # 计算各种匹配统计
    total = len(data)
    exact_matches = 0
    start_matches = 0
    end_matches = 0
    length_matches = 0
    no_matches = 0
    
    query_lengths = []
    ref_lengths = []
    
    for q_start, q_end, r_start, r_end in data:
        q_len = q_end - q_start
        r_len = r_end - r_start
        query_lengths.append(q_len)
        ref_lengths.append(r_len)
        
        if q_start == r_start and q_end == r_end:
            exact_matches += 1
        elif q_start == r_start:
            start_matches += 1
        elif q_end == r_end:
            end_matches += 1
        elif q_len == r_len:
            length_matches += 1
        else:
            no_matches += 1
    
    print(f"完全匹配 (位置+长度): {exact_matches:>3} ({exact_matches/total*100:.1f}%)")
    print(f"起始位置匹配:        {start_matches:>3} ({start_matches/total*100:.1f}%)")
    print(f"结束位置匹配:        {end_matches:>3} ({end_matches/total*100:.1f}%)")
    print(f"仅长度匹配:          {length_matches:>3} ({length_matches/total*100:.1f}%)")
    print(f"无匹配:              {no_matches:>3} ({no_matches/total*100:.1f}%)")
    
    # 长度分析
    avg_query_len = sum(query_lengths) / len(query_lengths)
    avg_ref_len = sum(ref_lengths) / len(ref_lengths)
    
    print(f"\n📏 长度分析:")
    print(f"平均Query长度: {avg_query_len:.2f}")
    print(f"平均Reference长度: {avg_ref_len:.2f}")
    print(f"长度差异: {abs(avg_query_len - avg_ref_len):.2f}")

def create_matplotlib_heatmap(data):
    """尝试创建matplotlib热力图"""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        query_starts = [item[0] for item in data]
        ref_starts = [item[2] for item in data]
        
        # 创建2D直方图
        plt.figure(figsize=(12, 10))
        
        # 使用hexbin创建热力图
        plt.hexbin(ref_starts, query_starts, gridsize=30, cmap='YlOrRd')
        plt.colorbar(label='匹配数量')
        
        plt.xlabel('Reference Start Position')
        plt.ylabel('Query Start Position')
        plt.title('Query Start vs Reference Start 热力图')
        
        # 添加对角线参考
        min_val = min(min(query_starts), min(ref_starts))
        max_val = max(max(query_starts), max(ref_starts))
        plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')
        plt.legend()
        
        plt.tight_layout()
        plt.show()
        
        return True
        
    except ImportError:
        print("matplotlib未安装，跳过图形化显示")
        return False
    except Exception as e:
        print(f"matplotlib显示失败: {e}")
        return False

# 主函数
if __name__ == "__main__":
    try:
        # 读取数据文件
        filename = 'data.txt'  # 可以修改为你的文件名
        print(f"正在读取文件: {filename}")
        data = parse_data_file(filename)
        print(f"成功读取 {len(data)} 条记录")
        
        # 创建热力图矩阵
        matrix, query_info, ref_info = create_heatmap_matrix(data, grid_size=30)
        
        # 打印ASCII热力图
        print_ascii_heatmap(matrix, query_info, ref_info)
        
        # 详细分析
        create_detailed_heatmap(data)
        
        # 模式分析
        analyze_matching_patterns(data)
        
        # 尝试matplotlib版本
        print(f"\n" + "="*60)
        print("尝试创建matplotlib热力图...")
        success = create_matplotlib_heatmap(data)
        
        if not success:
            print("如需图形化热力图，请安装matplotlib:")
            print("pip install matplotlib numpy")
        
    except FileNotFoundError:
        print(f"错误: 找不到文件 '{filename}'")
        print("请确保文件存在于当前目录中")
    except Exception as e:
        print(f"发生错误: {e}")
        import traceback
        traceback.print_exc()