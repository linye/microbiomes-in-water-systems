import os
from collections import defaultdict
import time

# 定义文件路径
cluster_file = "protein_clusters.txt"
group_file = "env.txt"
output_file = "cluster_env_sample_number_matrix.txt"

def process_group_file():
    """处理分组文件，返回样本到分组的映射和排序后的分组列表"""
    sample_to_group = {}
    groups = set()

    print("开始处理分组文件...")
    with open(group_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            sample, group = parts[0], parts[1]
            sample_to_group[sample] = group
            groups.add(group)

            if line_num % 1000000 == 0:
                print(f"已处理 {line_num} 行分组数据")

    sorted_groups = sorted(groups)
    print(f"分组处理完成，共发现 {len(sorted_groups)} 个分组")
    return sample_to_group, sorted_groups

def process_clusters(sample_to_group, groups):
    """处理聚类文件，正确提取代表序列"""
    print("开始处理聚类文件并生成结果...")
    start_time = time.time()
    processed_count = 0

    with open(output_file, 'w') as out_f:
        # 写入表头
        out_f.write("Representative\t" + "\t".join(groups) + "\n")

        with open(cluster_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                # 关键修改：正确分离代表序列和成员序列
                parts = line.split('\t')
                representative = parts[0]  # 总是第一列为代表序列
                
                # 提取所有序列（包括代表序列）
                all_sequences = [representative]
                if len(parts) > 1:
                    all_sequences.extend(parts[1].split(','))
                
                # 提取样本ID
                samples = set(seq.split('_')[0] for seq in all_sequences)
                
                # 统计每个分组的样本数
                group_counts = defaultdict(int)
                for sample in samples:
                    group = sample_to_group.get(sample)
                    if group:
                        group_counts[group] += 1
                
                # 准备输出行（使用正确的代表序列）
                counts = [str(group_counts.get(group, 0)) for group in groups]
                out_line = f"{representative}\t" + "\t".join(counts) + "\n"
                
                # 立即写入输出文件
                out_f.write(out_line)
                
                processed_count += 1
                if line_num % 100000 == 0:
                    elapsed = time.time() - start_time
                    rate = processed_count / elapsed if elapsed > 0 else 0
                    print(f"已处理 {line_num} 行, 速率: {rate:.2f} 行/秒, 耗时: {elapsed:.2f}秒")

    print(f"聚类处理完成，共处理 {processed_count} 个聚类")
    return processed_count

def main():
    # 步骤1: 处理分组文件
    sample_to_group, groups = process_group_file()

    # 步骤2: 处理聚类文件并直接输出结果
    cluster_count = process_clusters(sample_to_group, groups)

    print(f"结果已输出到文件 {output_file}")

if __name__ == "__main__":
    main()  # 修复：确保main()调用正确缩进
