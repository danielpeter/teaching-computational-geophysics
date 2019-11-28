--------------------------------
readme
--------------------------------

how to export a mesh for SPECFEM3D:

0. start Cubit/Trelis and create a meshed volume

1. run boundary script: menu bar 'Play Journal File'  -> select ./run_boundary_definition.py

2. export mesh: menu bar 'Play Journal File' -> select ./run_cubit2specfem.py

This will create all necessary files in folder MESH/

Move this folder to your SPECFEM example directory and setup the DATA/ folder with the corresponding 
Par_file, STATIONS and CMTSOLUTION files.


   
