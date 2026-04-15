import os
import numpy as np
from Bio.PDB import PDBParser

# 加载当前文件夹的所有PDB文件
def load_pdb_files():
    pdb_files = [f for f in os.listdir() if f.endswith('.pdb')]
    return pdb_files

# 计算两个PDB文件的RMSD（仅使用CA原子）
def calculate_rmsd(file1, file2):
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure('pdb1', file1)
    structure2 = parser.get_structure('pdb2', file2)

    coords1 = np.array([atom.get_coord() for atom in structure1.get_atoms() if atom.get_id() == 'CA'])
    coords2 = np.array([atom.get_coord() for atom in structure2.get_atoms() if atom.get_id() == 'CA'])

    if coords1.size == 0 or coords2.size == 0 or coords1.shape != coords2.shape:
        return np.inf

    diff = coords1 - coords2
    return np.sqrt(np.sum(diff ** 2) / len(coords1))

# 构建RMSD矩阵
def build_rmsd_matrix(pdb_files):
    num_files = len(pdb_files)
    rmsd_matrix = np.zeros((num_files, num_files))

    for i in range(num_files):
        for j in range(i + 1, num_files):
            rmsd = calculate_rmsd(pdb_files[i], pdb_files[j])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    # 用最大有效值替代 inf
    inf_mask = np.isinf(rmsd_matrix)
    if np.any(inf_mask):
        max_val = np.nanmax(rmsd_matrix[~inf_mask]) if np.any(~inf_mask) else 100.0
        rmsd_matrix[inf_mask] = max_val

    return rmsd_matrix

# 找到平均RMSD最小的文件
def find_center_by_avg_rmsd(rmsd_matrix, pdb_files):
    avg_distances = rmsd_matrix.mean(axis=1)
    center_index = np.argmin(avg_distances)
    return pdb_files[center_index]

# 重命名并只保留坐标部分的行
def rename_and_clean_pdb(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # 只保留以ATOM或HETATM开头的行
    filtered_lines = [line for line in lines if line.startswith('ATOM') or line.startswith('HETATM')]

    with open(output_file, 'w') as f:
        f.writelines(filtered_lines)

    os.remove(input_file)

# 主函数
def main():
    pdb_files = load_pdb_files()
    if len(pdb_files) == 0:
        print("No pdb file")
        return
    elif len(pdb_files) == 1:
        print("Only pdb file--> top1.pdb")
        rename_and_clean_pdb(pdb_files[0], "top1.pdb")
        return

    rmsd_matrix = build_rmsd_matrix(pdb_files)
    cluster_center = find_center_by_avg_rmsd(rmsd_matrix, pdb_files)
    rename_and_clean_pdb(cluster_center, "top1.pdb")
    print("Cluster center--> top1.pdb")

# 执行程序
if __name__ == "__main__":
    main()
