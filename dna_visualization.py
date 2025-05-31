import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import sys
from matplotlib.ticker import MultipleLocator
import edlib

def get_complement(base):
    """获取碱基的互补碱基"""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_dict.get(base, base)

def get_reverse_complement(seq):
    """Get reverse complement of a sequence"""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement_dict[base] for base in reversed(seq))

def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
    """Calculate edit distance between two sequences, considering both forward and reverse complement"""
    A = ref[ref_st:ref_en]
    a = query[query_st:query_en]
    _a = get_reverse_complement(query[query_st:query_en])
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

def find_matches(seq1, seq2, window_size=10):
    """
    查找两串DNA序列的匹配情况
    返回一个矩阵，其中：
    0: 不匹配
    1: 完全匹配
    2: 反向匹配
    3: 反向互补匹配
    4: 单核苷酸匹配
    5: 删除/插入
    """
    len1, len2 = len(seq1), len(seq2)
    match_matrix = np.zeros((len1, len2), dtype=int)
    
    # 完全匹配检测
    for i in range(len1 - window_size + 1):
        for j in range(len2 - window_size + 1):
            # 检查正向匹配
            if seq1[i:i+window_size] == seq2[j:j+window_size]:
                match_matrix[i:i+window_size, j:j+window_size] = 1
            
            # 检查反向匹配
            reverse_seq2 = seq2[j:j+window_size][::-1]
            if seq1[i:i+window_size] == reverse_seq2:
                match_matrix[i:i+window_size, j:j+window_size] = 2
            
            # 检查反向互补匹配
            reverse_complement_seq2 = get_reverse_complement(seq2[j:j+window_size])
            if seq1[i:i+window_size] == reverse_complement_seq2:
                match_matrix[i:i+window_size, j:j+window_size] = 3
    
    # 单核苷酸突变检测
    for i in range(len1):
        for j in range(len2):
            if match_matrix[i,j] == 0 and seq1[i] == seq2[j]:
                match_matrix[i,j] = 4
    
    # 删除/插入检测（通过查找较短的匹配片段）
    for i in range(len1 - 5 + 1):
        for j in range(len2 - 5 + 1):
            if match_matrix[i:i+5, j:j+5].sum() >= 3:  # 如果5个碱基中有3个以上匹配
                match_matrix[i:i+5, j:j+5] = np.where(match_matrix[i:i+5, j:j+5] == 0, 5, match_matrix[i:i+5, j:j+5])
    return match_matrix

def plot_dna_alignment(seq1, seq2, match_matrix, additional_regions=None):
    """
    绘制DNA序列匹配情况的可视化图
    additional_regions: 列表，每个元素为 (query_start, query_end, length, ref_start, ref_end, reversed) 的元组
    """
    # 设置更大的图形尺寸和更高的DPI
    plt.figure(figsize=(20, 20), dpi=300)
    
    # 创建自定义颜色映射 - 使用更鲜明的颜色
    colors = ['#FFFFFF',  # 白色 - 不匹配
              '#FF0000',  # 红色 - 完全匹配
              '#0000FF',  # 蓝色 - 反向匹配
              '#800080',  # 紫色 - 反向互补匹配
              '#00FF00',  # 绿色 - 单核苷酸匹配
              '#FFD700']  # 金色 - 删除/插入
    cmap = ListedColormap(colors)
    
    # 创建主图
    ax = plt.gca()
    im = ax.imshow(match_matrix, cmap=cmap, aspect='auto', interpolation='nearest')
    
    # 添加网格线
    ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # 设置刻度
    major_interval = 100  # 主刻度间隔
    minor_interval = 20   # 次刻度间隔
    
    # 设置主刻度和次刻度
    ax.xaxis.set_major_locator(MultipleLocator(major_interval))
    ax.xaxis.set_minor_locator(MultipleLocator(minor_interval))
    ax.yaxis.set_major_locator(MultipleLocator(major_interval))
    ax.yaxis.set_minor_locator(MultipleLocator(minor_interval))
    
    # 添加额外的区域标记
    if additional_regions:
        for region in additional_regions:
            query_st, query_en, length, ref_st, ref_en, is_reversed = region
            # 使用红色斜线标记这些区域
            if is_reversed:
                # 如果是反向的，画反向的斜线
                ax.plot([query_st, query_en], [ref_en, ref_st], 
                       color='#FF4500', linewidth=5, linestyle='--', alpha=0.7)
            else:
                # 如果是正向的，画正向的斜线
                ax.plot([query_st, query_en], [ref_st, ref_en], 
                       color='#FF4500', linewidth=5, linestyle='--', alpha=0.7)
    
    # 设置坐标轴标签
    plt.xlabel('Sequence 2 Position', fontsize=12, fontweight='bold')
    plt.ylabel('Sequence 1 Position', fontsize=12, fontweight='bold')
    
    # 添加图例
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor='#FFFFFF', edgecolor='black', label='No Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FF0000', label='Perfect Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#0000FF', label='Reverse Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#800080', label='Reverse Complement Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#00FF00', label='Single Nucleotide Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FFD700', label='Indel/Short Match')
    ]
    
    # 如果有额外区域，添加图例
    if additional_regions:
        legend_elements.append(plt.Line2D([0], [0], color='#FF4500', linestyle='--', 
                                        label='Additional Regions', linewidth=2))
    
    # 将图例放在图形外部
    plt.legend(handles=legend_elements, 
              loc='center left', 
              bbox_to_anchor=(1.05, 0.5),
              fontsize=10,
              frameon=True,
              edgecolor='black',
              title='Match Types',
              title_fontsize=12)
    
    # 设置标题
    plt.title('DNA Sequence Alignment Visualization\n' + 
              f'Sequence 1 Length: {len(seq1)}, Sequence 2 Length: {len(seq2)}',
              fontsize=14, 
              pad=20)
    
    # 调整布局以确保图例完全显示
    plt.tight_layout()
    
    # 保存图片，确保图例不被裁剪
    plt.savefig('dna_alignment.png', 
                dpi=300, 
                bbox_inches='tight',
                facecolor='white',
                edgecolor='none')
    plt.close()

def find_matching_regions(match_matrix, min_length=10):
    """
    从匹配矩阵中找出所有匹配区间
    返回格式为 (query_start, query_end, ref_start, ref_end, match_type) 的列表
    match_type: 1=完全匹配, 2=反向匹配, 3=反向互补匹配
    """
    regions = []
    rows, cols = match_matrix.shape
    
    # 对每种匹配类型分别处理
    for match_type in [1, 2, 3]:
        # 创建当前匹配类型的掩码
        mask = (match_matrix == match_type)
        
        # 在行方向上查找连续匹配
        for i in range(rows):
            start = None
            for j in range(cols):
                if mask[i, j] and start is None:
                    start = j
                elif not mask[i, j] and start is not None:
                    if j - start >= min_length:
                        regions.append((i, i, start, j-1, match_type))
                    start = None
            if start is not None and cols - start >= min_length:
                regions.append((i, i, start, cols-1, match_type))
        
        # 在列方向上查找连续匹配
        for j in range(cols):
            start = None
            for i in range(rows):
                if mask[i, j] and start is None:
                    start = i
                elif not mask[i, j] and start is not None:
                    if i - start >= min_length:
                        regions.append((start, i-1, j, j, match_type))
                    start = None
            if start is not None and rows - start >= min_length:
                regions.append((start, rows-1, j, j, match_type))
    
    # 按匹配类型和起始位置排序
    regions.sort(key=lambda x: (x[4], x[0], x[2]))
    return regions

def find_optimal_regions(match_matrix, seq1, seq2, min_length=30, max_edit_ratio=0.1):
    """
    Find optimal matching regions that maximize the score
    Returns a list of (query_start, query_end, ref_start, ref_end) tuples
    """
    regions = []
    rows, cols = match_matrix.shape
    used_positions = set()  # Track used positions to avoid overlaps
    
    # Process each type of match (1, 2, 3) separately
    for match_type in [1, 2, 3]:
        mask = (match_matrix == match_type)
        
        # Find continuous regions
        for i in range(rows):
            for j in range(cols):
                if mask[i, j] and (i, j) not in used_positions:
                    # Try to extend the region
                    length = 1
                    while (i + length < rows and j + length < cols and 
                           mask[i + length, j + length] and 
                           (i + length, j + length) not in used_positions):
                        length += 1
                    
                    if length >= min_length:
                        # Calculate edit distance
                        edit_dist = calculate_distance(seq1, seq2, j, j+length, i, i+length)
                        if edit_dist / length <= max_edit_ratio:
                            # Mark positions as used
                            for k in range(length):
                                used_positions.add((i + k, j + k))
                            # Only store the four coordinates
                            regions.append((i, i+length-1, j, j+length-1))
    
    # Sort regions by length (descending) to prioritize longer matches
    regions.sort(key=lambda x: (x[1]-x[0]), reverse=True)
    return regions

def calculate_score(regions, seq1, seq2):
    """Calculate the total score for the given regions"""
    total_score = 0
    for region in regions:
        query_st, query_en, ref_st, ref_en = region
        length = query_en - query_st + 1
        edit_dist = calculate_distance(seq1, seq2, ref_st, ref_en+1, query_st, query_en+1)
        total_score += max(length - edit_dist, 0)
    return total_score

def print_matching_regions(regions, seq1, seq2):
    """Print matching regions and their scores"""
    print("\nMatching Regions (query_start, query_end, ref_start, ref_end):")
    print("-" * 80)
    total_score = 0
    for region in regions:
        query_st, query_en, ref_st, ref_en = region
        length = query_en - query_st + 1
        edit_dist = calculate_distance(seq1, seq2, ref_st, ref_en+1, query_st, query_en+1)
        score = max(length - edit_dist, 0)
        total_score += score
        print(f"({query_st}, {query_en}, {ref_st}, {ref_en}) - Length: {length}, Edit Distance: {edit_dist}, Score: {score}")
    print(f"\nTotal Score: {total_score}")

def main():
    # Read input sequences
    print("Please enter the first DNA sequence:")
    seq1 = sys.stdin.readline().strip().upper()
    print("Please enter the second DNA sequence:")
    seq2 = sys.stdin.readline().strip().upper()
    print(seq1)
    print(seq2)
    # Validate input
    valid_bases = {'A', 'T', 'C', 'G'}
    if not all(base in valid_bases for base in seq1) or not all(base in valid_bases for base in seq2):
        print("Error: Invalid bases in input sequence. Only A, T, C, G are allowed.")
        return
    
    # Read additional regions
    print("\nEnter additional regions (one per line, format: query_start query_end length ref_start ref_end reversed)")
    print("Enter an empty line to finish input:")
    additional_regions = []
    while True:
        line = sys.stdin.readline().strip()
        print(line)
        if not line:
            break
        try:
            # Parse the input line
            parts = line.split()
            if len(parts) != 6:
                print(f"Error: Invalid input format. Skipping line: {line}")
                continue
            query_st, query_en, length, ref_st, ref_en, reversed_flag = map(int, parts)
            additional_regions.append((query_st, query_en, length, ref_st, ref_en, bool(reversed_flag)))
        except ValueError as e:
            print(f"Error parsing line: {line}. Skipping...")
            continue
    print(len(additional_regions))
    print("Next ")

    # Find matches
    match_matrix = find_matches(seq1, seq2)
    print("Next1")
    # Find optimal regions
    regions = find_optimal_regions(match_matrix, seq1, seq2)
    print("Next2")
    # Print regions and scores
    print_matching_regions(regions, seq1, seq2)
    print("Next3")
    # Plot visualization with additional regions
    plot_dna_alignment(seq1, seq2, match_matrix, additional_regions)
    print("\nVisualization has been saved as 'dna_alignment.png'")

if __name__ == "__main__":
    main() 