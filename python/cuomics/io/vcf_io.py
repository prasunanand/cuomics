import vcf
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

def geno_read(file):
  vcf_reader = vcf.Reader(open(file, 'r'))
  df = pd.DataFrame(columns = ['Chrom', 'Pos', 'Ref', 'Alt', 'GT'])
  for record in vcf_reader:
    gt_samples = []
    for call in record:
      sample_id = call.sample
      gt = convert_gt_to_score(call['GT'])
      gt_samples.append(gt)
    df = df.append(pd.DataFrame([[str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT), gt_samples]], columns=['Chrom', 'Pos', 'Ref', 'Alt', 'GT']))
  return df

def convert_gt_to_score(gt):
  # print(j)
  jk = gt.split("|")
  if len(jk) != 2:
    jk = gt.split("/")
  if(jk[0] == '.'):
    return 0
  elif(jk[1] == '.'):
    return 0
  return get_gt(int(jk[0]),int(jk[1]))

def get_gt(j, k):
  return k*(k+1)/2 + j

def store_as_parquet(df):
  table = pa.Table.from_pandas(df)
  pq.write_table(table, 'snps.parquet')

geno = geno_read('./python/cuomics/tests/example.vcf')
print(geno)
store_as_parquet(geno)
