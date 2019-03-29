import numpy as np
from cuomics.io.parquet_io import read_from_parquet

def get_genotype(raw_genotype_matrix, number_of_samples, biallelic_autosomal_variants):
  genotype_matrix = np.empty([number_of_samples, biallelic_autosomal_variants])
  for i in range(number_of_samples):
    for j in range(biallelic_autosomal_variants):
      p = np.mean(genotype_matrix[:, j])/2
      genotype_matrix[i,j] = (raw_genotype_matrix[i,j] - 2*p)/np.sqrt(2*p*(1-p)*biallelic_autosomal_variants)

  return genotype_matrix

def get_kinship(genotype_matrix):
  # Kinship/Correlation/Covariance Matrix
  kinship_matrix = np.dot(genotype_matrix, genotype_matrix.T)
  return kinship_matrix

def pca(kinship_matrix):
  # kinship_matrix=U*S*V.T
  U, S, V = np.linalg.svd(kinship_matrix)
  print(U)
  return U,S,V

def pca(kinship_matrix):
  # M=U*S*V.T
  # columns of U are left singular vectors (orthonormal in Rn),
  # columns of V are right singular vectors (orthonormal in Rm),
  # and S=diag(s1,s2,…) with ordered singular values s1≥s2≥⋯≥0.
  U, S, VT = np.linalg.svd(kinship_matrix)
  return U, S, VT

def main():
  df = read_from_parquet('../snps.parquet')

  number_of_samples = 3
  biallelic_autosomal_variants = 6

  raw_genotype_matrix  = np.empty([number_of_samples, biallelic_autosomal_variants])
  i = 0


  for row in df.GT:
    raw_genotype_matrix[:, i] = np.array(row)
    i+=1

  print("============raw_genotype_matrix==============")
  print(raw_genotype_matrix)

  print("============genotype_matrix==================")
  genotype_matrix = get_genotype(raw_genotype_matrix, number_of_samples, biallelic_autosomal_variants)
  print(genotype_matrix)

  print("============kinship_matrix===================")
  kinship_matrix = get_kinship(genotype_matrix)
  print(kinship_matrix)

  print("============PCA_analysis=====================")
  U, S, VT = pca(kinship_matrix)
  print(U)
  print(S)
  print(VT)

if __name__ == '__main__':
  main()