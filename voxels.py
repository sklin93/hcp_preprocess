import nibabel as nib
import numpy as np
mask_file = '/share/igert/share/DSI_data/voxel_clustering_MDA/MNI_mask_MC.nii.gz'
mask = nib.load(mask_file).get_data() > 0

idx2voxelIJK = np.array(np.where(mask)).T

print 'Dimension:', idx2voxelIJK.shape, '\n'
print idx2voxelIJK
print idx2voxelIJK[0]
