# crystal
A simple C++ program to create arbitrarily shaped, crystalline atomic clusters. 

This program takes a crystal structure definition (unit cell and basis), creates an na x nb x nc supercell, and fills it with atoms according to the crystal structure. Then it carves out a cluster by removing all atoms that fall out of the volume bounded by a set of planes specified in the input. Out is defined as the possitive side of the normal vector defining each plane.

# Usage
```
crystal.exe input_file.in output_file.xyz
```
# Input

A text file with the following structure:
* ```Name of structure to be created```a text string that will be written in the output.
* The text string: ```Scale factor in angsroms```
* A decimal number, setting the unit length of the basis vectors
* A blank line
* The text string: ```Size```
* ```na nb  nc``` Three tab-separated integers specifying the size of the supercell to be created (number of unit cells along each unit vector).
* A blank line
* The text string: ```Direct lattice vectors (Cartesian canonic coordinates)```
* Three lines, each contaiing three tab-separated decimal numbers specifying the coordinates of a unit cell basis vector. 
* A blank line
* The text string: ```Basis vectors```
* ```n_basis``` An integer specifying the number of atoms in the crystal basis.
* n_basis lines with the format ```Sym  x1  x2  x3``` (tab-separated); where Sym is the atomic symbol of the atom sitting at (reduced coordinates) (x1, x2, x3).   
* A blank line
* The text string: ```Planes```
* ```n_planes``` An integer specifying the number of planes to be carved out of the cluster
* n_planes groups of two lines, each containing three tab_separated numbers:
  * plane normal vector
  * Coordinates of a point in the plane, to set distance from origin


# Output
An [xyz file](http://wiki.jmol.org/index.php/File_formats/Formats/XYZ) containing the coordinates of all atoms in the cluster, readable by visulaiztion programs; e.g. [Jmol](http://jmol.sourceforge.net). The descriptive string in the xyz file is the ```Name of structure to be created``` string in the input file.

# Further documentation
\exmples\ contains some examples of input (.in) and output (.xyz) files. Each input and output file name differs only in the extension.
\images\ conatins some Jmol-generated images of clusters created with crytal 
