import pandas as pd

df_base = './mda_fibs/%s/warped_mda.h5'
load_df = lambda subj: pd.read_hdf(df_base % subj, str(subj))

subj = '180129'

scan = load_df(subj)
print('Number of voxels:', len(scan))
print(scan.head())

scan.to_csv('warped_mda_180129.csv', sep='\t', index=False)

