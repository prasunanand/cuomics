import pyarrow as pa
import pandas as pd
import pyarrow.parquet as pq

def read_from_parquet(file_name = 'snps.parquet'):
  table2 = pq.read_table(file_name)
  return table2.to_pandas()

# geno = read_from_parquet('snps.parquet')
# print(geno)
# geno.GT = geno.GT + 1
# print(geno.GT + 1)