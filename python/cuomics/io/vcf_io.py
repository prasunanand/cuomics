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
      # print(sample_id, " -> ", gt)
      gt_samples.append(gt)
    df = df.append(pd.DataFrame([[str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT), gt_samples]], columns=['Chrom', 'Pos', 'Ref', 'Alt', 'GT']))
  return df

def convert_gt_to_score(gt):
  if(gt ==  '0|0' or gt ==  '0/0'):
    return 0
  elif(gt ==  '0|1' or gt ==  '0/1' or gt ==  '1|0' or gt ==  '1/0'):
    return 0
  elif(gt ==  '1|1' or gt ==  '1/1'):
    return 0
  else:
    # check how to handle
    # return 'NA'
    return 0

def store_as_parquet(df):
  table = pa.Table.from_pandas(df)
  pq.write_table(table, 'snps.parquet')

geno = geno_read('./python/cuomics/tests/example.vcf')
print(geno)
store_as_parquet(geno)
