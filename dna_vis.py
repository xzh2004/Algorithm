import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import sys

def find_matches(seq1, seq2, window_size=10):
    """
    查找两串DNA序列的匹配情况
    返回一个矩阵，其中：
    0: 不匹配
    1: 完全匹配
    2: 倒位匹配
    3: 单核苷酸突变
    4: 删除/插入
    """
    len1, len2 = len(seq1), len(seq2)
    match_matrix = np.zeros((len1, len2), dtype=int)
    
    # 完全匹配检测
    for i in range(len1 - window_size + 1):
        for j in range(len2 - window_size + 1):
            # 检查正向匹配
            if seq1[i:i+window_size] == seq2[j:j+window_size]:
                match_matrix[i:i+window_size, j:j+window_size] = 1
            
            # 检查倒位匹配
            reverse_seq2 = seq2[j:j+window_size][::-1]
            if seq1[i:i+window_size] == reverse_seq2:
                match_matrix[i:i+window_size, j:j+window_size] = 2
    
    # 单核苷酸突变检测
    for i in range(len1):
        for j in range(len2):
            if match_matrix[i,j] == 0 and seq1[i] == seq2[j]:
                match_matrix[i,j] = 3
    
    # 删除/插入检测（通过查找较短的匹配片段）
    for i in range(len1 - 5 + 1):
        for j in range(len2 - 5 + 1):
            if match_matrix[i:i+5, j:j+5].sum() >= 3:  # 如果5个碱基中有3个以上匹配
                match_matrix[i:i+5, j:j+5] = np.where(match_matrix[i:i+5, j:j+5] == 0, 4, match_matrix[i:i+5, j:j+5])
    
    return match_matrix

def plot_dna_alignment(seq1, seq2, match_matrix):
    """
    绘制DNA序列匹配情况的可视化图
    """
    # 创建自定义颜色映射
    colors = ['white', 'red', 'blue', 'green', 'yellow']
    cmap = ListedColormap(colors)
    
    # 创建图形
    plt.figure(figsize=(12, 12))
    plt.imshow(match_matrix, cmap=cmap, aspect='auto')
    
    # 设置坐标轴标签
    plt.xlabel('Sequence 2 Position')
    plt.ylabel('Sequence 1 Position')
    
    # 添加图例
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', label='No Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='red', label='Perfect Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='blue', label='Inversion'),
        plt.Rectangle((0, 0), 1, 1, facecolor='green', label='Single Nucleotide Match'),
        plt.Rectangle((0, 0), 1, 1, facecolor='yellow', label='Indel/Short Match')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    # 设置标题
    plt.title('DNA Sequence Alignment Visualization')
    
    # 保存图片
    plt.savefig('dna_alignment.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # 读取输入序列
    print("请输入第一串DNA序列：")
    seq1 = sys.stdin.readline().strip().upper()
    print("请输入第二串DNA序列：")
    seq2 = sys.stdin.readline().strip().upper()
    
    # 验证输入
    valid_bases = {'A', 'T', 'C', 'G'}
    if not all(base in valid_bases for base in seq1) or not all(base in valid_bases for base in seq2):
        print("错误：输入序列包含无效的碱基。只允许 A, T, C, G。")
        return
    
    # 查找匹配
    match_matrix = find_matches(seq1, seq2)
    
    # 绘制图形
    plot_dna_alignment(seq1, seq2, match_matrix)
    print("可视化结果已保存为 'dna_alignment.png'")

if __name__ == "__main__":
    main() 