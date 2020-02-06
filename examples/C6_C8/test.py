import mdtraj as md
import auxerre as axr
import numpy as np

pdb_c6 = md.load('C6.pdb')
pdb_c8 = md.load('C8.pdb')

traj_c6 = md.iterload('1.xyz.gz', top=pdb_c6, chunk=1000)
traj_c8 = md.iterload('2.xyz.gz', top=pdb_c8, chunk=1000)

analyzer = axr.Analyzer(
    axr.DiffusivityCalculator(),
    pdb_c6.unitcell_lengths,
)

for i, (chunk_c6, chunk_c8) in enumerate(zip(traj_c6, traj_c8)):
    chunk = chunk_c8.stack(chunk_c6)
    print(i, chunk)
    analyzer.add_trajectory(chunk)
    # if i == 0:
    #     break
analyzer.correlate()
