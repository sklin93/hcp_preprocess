# hcp_preprocess
Original codes credit to Tanyi.

Changed:
[vector_bin.py](vector_bin.py) for mapping whitematter voxels into n*n*n blocks, which makes later computation can be fed into memory.
[ER.py](ER.py) for calculating structural matrix of grey matter regions based on effective resistance.

Running order:
wm.py & gm.py
vector_bin.py
ER.py
