import numpy as np
from sklearn.cluster import KMeans
from Bio import PDB

# Function to read multiple conformations from a PDB file
def read_conformations(pdb_file):
    conformations = []
    parser = PDB.PDBParser(QUIET=True)

    with open(pdb_file, 'r') as file:
        content = file.read()
        structures = content.split('END')  # Split structures by the "END" keyword

    for idx, structure_str in enumerate(structures):
        if structure_str.strip():  # Ignore empty splits
            filename = f"temp_structure_{idx}.pdb"
            with open(filename, 'w') as temp_file:
                temp_file.write(structure_str)

            try:
                structure = parser.get_structure(f"conf_{idx}", filename)
                conformations.append(structure)
            except Exception as e:
                print(f"Error parsing conformation {idx}: {e}")

    return conformations

# Function to calculate RMSD between two structures
def calculate_rmsd(structure1, structure2):
    coords1 = np.array([atom.get_coord() for atom in structure1.get_atoms() if atom.get_id() == "P"])
    coords2 = np.array([atom.get_coord() for atom in structure2.get_atoms() if atom.get_id() == "P"])

    if len(coords1) != len(coords2):
        raise ValueError("Structures have differing numbers of atoms")

    diff = coords1 - coords2
    return np.sqrt(np.sum(diff ** 2) / len(coords1))

# Generate feature matrix for clustering
def generate_feature_matrix(conformations):
    features = []
    for conformation in conformations:
        coords = [atom.get_coord() for atom in conformation.get_atoms() if atom.get_id() == "P"]
        features.append(np.array(coords).flatten())
    return np.array(features)

# Main clustering function
def cluster_conformations(pdb_file, n_clusters=3):
    conformations = read_conformations(pdb_file)
    feature_matrix = generate_feature_matrix(conformations)

    # Apply clustering algorithm (e.g., K-means)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    labels = kmeans.fit_predict(feature_matrix)

    # Output results
    for idx, label in enumerate(labels):
        print(f"Conformation {idx}: Cluster {label}")

    return labels

# Run the clustering process
if __name__ == "__main__":
    pdb_file = "cf.pdb"  # Path to your PDB file
    n_clusters = 5       # Number of clusters to form

    labels = cluster_conformations(pdb_file, n_clusters)

    with open("clustering_results.txt", "w") as result_file:
        for idx, label in enumerate(labels):
            result_file.write(f"Conformation {idx}: Cluster {label}\n")

    print("Clustering complete. Results saved to clustering_results.txt")

