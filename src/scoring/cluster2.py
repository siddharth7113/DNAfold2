import os
import shutil

# 文件和文件夹路径
clustering_file = "clustering_results.txt"
source_folder = "."

def main():
    # 读取 clustering_results.txt 内容
    with open(clustering_file, "r") as file:
        lines = file.readlines()

    # 解析 cluster 信息
    cluster_map = {}
    for line in lines:
        parts = line.strip().split(":")
        if len(parts) == 2:
            conformation = parts[0].strip().split()[1]  # 获取 Conformation n 的 n
            cluster = parts[1].strip().split()[1]      # 获取 Cluster 的编号
            cluster_map[conformation] = cluster

    # 创建 cluster 文件夹并移动 .pdb 文件
    for conformation, cluster in cluster_map.items():
        source_file = os.path.join(source_folder, f"temp_structure_{conformation}.pdb")
        target_folder = os.path.join(source_folder, f"Cluster_{cluster}")

        # 创建目标文件夹（如果不存在）
        os.makedirs(target_folder, exist_ok=True)

        # 移动文件
        if os.path.exists(source_file):
            shutil.move(source_file, os.path.join(target_folder, f"temp_structure_{conformation}.pdb"))
            print(f"Moved {source_file} to {target_folder}")
        else:
            print(f"File {source_file} not found. Skipping.")

if __name__ == "__main__":
    main()

