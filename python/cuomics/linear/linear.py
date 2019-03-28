from ..io import read_geno, read_pheno, read_covariates, read_variants, read_annotations
from ..preprocessing import process_covariates_phenotypes


def lm_run(variants_file, pheno_file, covariates_file, geno_file, bfile, mode):
  print("Running Linear Regression")

  phenotypes, indicator_phenotypes  = read_pheno(pheno_file)

  covariates, indicator_covariates, num_covariates = read_covariates(covar_file)

  num_test_individuals, indicator_individuals = process_covariates_phenotypes(indicator_phenotypes,
                                                                              covariates,
                                                                              indicator_covariates,
                                                                              num_covariates)

  num_total_individuals = len(indicator_individuals)

  # rename later to set_variants
  set_variants = read_variants(variants_file)

  # mappers/dictionary

  mapRS2chr, mapRS2bp, mapRS2cm = read_annotations(variants_file)

  num_phenotypes = 1

  Y = np.zeros(num_test_individuals, num_phenotypes)
  W = np.zeros(num_test_individuals, num_covariates)

  copy_covariates_phenotypes(W, Y,
                            indicator_individuals, indicator_covariates,
                            phenotypes, covariates,
                            num_phenotypes, num_covariates)

  indicator_variants, variant_info = read_geno(geno_file, num_total_individuals, W,
                                       indicator_individuals, set_variants,
                                       mapRS2cm, mapRS2bp, mapRS2chr)

  WtW = np.dot(W.T, W)
  WtWi = np.inverse(WtW)

  WtY = np.dot(W.T, Y)

  summary_of_statistics = analyze_bimbam(geno_file, W, Y, WtWi, WtY,
                                         indicator_individuals, indicator_variants,
                                         num_test_individuals, num_total_individuals);

  return summary_of_statistics
