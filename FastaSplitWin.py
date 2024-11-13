from Bio import SeqIO
import argparse

# 创建命令行参数解析器
parser = argparse.ArgumentParser(description="FastaSplitWin -将fasta文件分割成多个窗口，其结果每个窗口文件都包含原文件的所有序列。")
parser.add_argument("-i","--input_fasta", type=str,required=True,help="输入fasta文件路径")
parser.add_argument("-w","--window_size", type=int, default=5000, help="分割窗口大小 (default: 5000)")
parser.add_argument("-o","--output", type=str, default="output", help="输出文件路径（文件名），不用额外添加后缀。 (default: 'output')")

# 解析命令行参数
args = parser.parse_args()

# 获取输入的FASTA文件路径
input_fasta = args.input_fasta

# 设置窗口大小
window_size = args.window_size

# 读取并分割FASTA文件
def split_fasta_by_window(input_fasta, window_size):
    # 读取输入的FASTA文件
    records = list(SeqIO.parse(input_fasta, "fasta"))
    seq_lengths = [len(record.seq) for record in records]
    
    # 计算最大窗口数量（基于最长的序列）
    max_windows = max((seq_len + window_size - 1) // window_size for seq_len in seq_lengths)

    # 使用一个字典来存储所有分割的窗口数据
    window_buffers = {f"{args.output}_window_{i+1}.fasta": [] for i in range(max_windows)}

    # 遍历每个序列并填充到对应窗口的buffer中
    for record in records:
        sequence = str(record.seq)
        seq_len = len(sequence)
        num_windows = (seq_len + window_size - 1) // window_size  # 当前序列的窗口数量

        for window_num in range(num_windows):
            start = window_num * window_size
            end = min(start + window_size, seq_len)  # 防止超出序列长度
            split_seq = sequence[start:end]
            window_id = f"{args.output}_window_{window_num + 1}.fasta"
            window_buffers[window_id].append(f">{record.id}_window_{window_num + 1}\n{split_seq}\n")

    # 将缓冲区内容写入对应的文件
    for window_file, content in window_buffers.items():
        with open(window_file, "w") as out_f:
            out_f.writelines(content)

# 运行分割函数
split_fasta_by_window(input_fasta, window_size)
