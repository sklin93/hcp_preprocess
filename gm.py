import numpy as np
import nibabel as nib
atlas_loc = './ROI_scale60.nii.gz'
atlas = nib.load(atlas_loc).get_data()

print('Atlas dimensions:', atlas.shape, '\n')
print('Unique region numbers:\n', np.unique(atlas), '\n\n')
print(np.count_nonzero(atlas))
# print(atlas)
atlas.tofile("atlas_scale60.csv", sep=',', format='%10.5f')
