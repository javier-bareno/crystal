//---------------------------------------------------------------------------

#pragma hdrstop

//---------------------------------------------------------------------------

#pragma argsused
#include <iostream>
#include <fstream>
#include <string>
#include <math>

using namespace std;

struct atom {
	string element;
	double x,y,z;
}; // this structure defines an atom as an element at some coordinates

struct many_atoms {
	int number; // number of atoms
	string name; // name of the crystal
	atom *atoms; // array of atoms, not allocated
};

struct plane { // defined by normal vector and origin
	double nx, ny, nz;
	double ox, oy, oz;
};

struct many_planes {
	int number; // number of planes
	plane *planes; // array to store plains, not allocated yet
};

struct crys_system {
	many_atoms *crystal;
	many_planes *cuts;
};


crys_system readdata(char* fname)
{
	//reads data from file and constructs 'cube' of atoms
	// version using striuctures

	double scale; // length of vector (1 0 0) in angstroms
	int maxi, maxj, maxk; // define size of initial crystal

	// direct lattice vectors' components
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	int nb; // number of basis vectors
	atom *basis; // array of atoms to store basis, not allocated yet
	atom *crystal; // array of atoms to store crystal, not allocated yet
	plane *planes; // array of planes to cut crystal
	string name; // name of structure, to be put on out file
	string istr; // input string
	double rx, ry, rz; // crystal point position, will be swept to create cryst
	int nplanes; // number of cutting planes
	int natoms, cnt;

	crys_system answer; // returns crystal and planes

	ifstream fin(fname);
	ofstream fout;

	getline(fin, name);// reads name of structure
	getline(fin, istr);// reads label
	fin >> scale;
	getline(fin, istr);// reads end of line
	getline(fin, istr);// reads blank line
	getline(fin, istr);// reads label

	fin >> maxi >>  maxj >>  maxk ;
	getline(fin, istr);// reads end of line
	getline(fin, istr);// reads blank line
	getline(fin, istr);// reads label

	fin >> ax >> ay >> az;
	getline(fin, istr);// reads end of line
	fin >> bx >> by >> bz;
	getline(fin, istr);// reads end of line
	fin >> cx >> cy >> cz;
	getline(fin, istr);// reads end of line
	getline(fin, istr);// reads blank line
	getline(fin, istr);// reads label

	fin >> nb;
	getline(fin, istr);// reads end of line

	// allocate memory for the basis of the crystal
	basis = new atom[nb];

	for (int i=0; i< nb; i++) {
		// read each basis vector
		fin >> basis[i].element >> basis[i].x >> basis[i].y >> basis[i].z;
		getline(fin, istr);// reads end of line
		// now scale atom coordinates
		basis[i].x *= scale;
		basis[i].y *= scale;
		basis[i].z *= scale;
	};
	getline(fin, istr);// reads blank line

	getline(fin, istr);// reads label 'planes'
	fin >> nplanes ; // reads number of planes
	planes = new plane[nplanes];
	getline(fin, istr);// reads end of line
	cout << istr;

	for (int i =0; i< nplanes; i++) {
		fin >> planes[i].nx;
		fin >> planes[i].ny;
		fin >> planes[i].nz;
		getline(fin, istr);// reads end of line
		fin >> planes[i].ox; planes[i].ox *= scale;
		fin >> planes[i].oy; planes[i].oy *= scale;
		fin >> planes[i].oz; planes[i].oz *= scale;
		getline(fin, istr);// reads end of line
	};
	fin.close();

	// generate crystal
	natoms = (maxi*maxj*maxk*nb);
	crystal = new atom[natoms]; // allocate memory to store crystal
	cnt = 0;
	for (int i=0; i<maxi; i++)
	{	for (int j=0; j<maxj; j++)
		{	for (int k=0; k<maxk; k++)
			{
				rx = scale * (i*ax + j*bx + k*cx);
				ry = scale * (i*ay + j*by + k*cy);
				rz = scale * (i*az + j*bz + k*cz);
				for (int l=0; l<nb; l++)
				{
					crystal[cnt].element= basis[l].element;
					crystal[cnt].x = basis[l].x + rx;
					crystal[cnt].y = basis[l].y + ry;
					crystal[cnt].z = basis[l].z + rz;
					cnt ++;
				};
			};
		};
	};
	//free memory storing basis
	delete[] basis;

	// prepare answer
	answer.crystal = new many_atoms;
	answer.cuts = new many_planes;
	answer.crystal->number=natoms;
	answer.cuts->number=nplanes;
	answer.crystal->name=name;
	answer.crystal->atoms=crystal;
	answer.cuts->planes=planes;

	return(answer);
};

void tofile(char* fname2, many_atoms *crystal)
{
	// copy crystal to out file
	ofstream fout(fname2);

	fout << crystal->number << "\n";  // writes size of structure
	fout << crystal->name << "\n"; // writes name of structure

	for (int cnt = 0; cnt < crystal->number; cnt++)
	{
		fout << crystal->atoms[cnt].element << "\t";
		fout << crystal->atoms[cnt].x  << "\t";
		fout << crystal->atoms[cnt].y << "\t";
		fout << crystal->atoms[cnt].z << "\n";
	};

	fout.close();

	return;
};

double distance(atom one_atom, plane one_plane)
{
	double n, nx, ny, nz;
	double no, ox, oy, oz;
	double nr, rx, ry, rz;

	nx = one_plane.nx;
	ny = one_plane.ny;
	nz = one_plane.nz;
	ox = one_plane.ox;
	oy = one_plane.oy;
	oz = one_plane.oz;
	rx = one_atom.x;
	ry = one_atom.y;
	rz = one_atom.z;
	n = sqrt(nx*nx + ny*ny + nz*nz);
	no = nx*ox + ny*oy + nz*oz;
	nr = nx*rx + ny*ry + nz*rz;

	return((nr-no)/n);
};


many_atoms slash(many_atoms crystal, plane cut)
{
	// cuts the crystal with a plane
	// returns cut crystal

	many_atoms cut_crystal;
	atom *atoms = new atom[crystal.number];
	int natoms = 0; // counts atoms remaining after cut

	for (int i=0; i<crystal.number; i++)
	{	if (distance(crystal.atoms[i], cut) < 0.1)
		{	atoms[natoms]=crystal.atoms[i];
			natoms++;
		};
	};
	cut_crystal.number = natoms;
	cut_crystal.name = crystal.name;
	if (natoms > 0)
	{	cut_crystal.atoms = new atom[natoms];
		for (int i=0; i<natoms; i++)
			cut_crystal.atoms[i] = atoms[i];
	};
	delete[] atoms;
	return (cut_crystal);
};

int main(int argc, char* argv[])
{
	// argv[1] points to input file
	// argv[2] points to output file

	many_atoms *crystal_cell;
	many_atoms crystal_to_cut, crystal_cut;
	many_planes *crystal_cuts;
	crys_system crystal;

	crystal = readdata(argv[1]); //creates crystal
	crystal_cell = crystal.crystal;
	crystal_cuts = crystal.cuts;

	if (crystal_cuts->number < 1)
	{	// no cuts
		tofile(argv[2], crystal_cell);
		//free memory
		 delete[] crystal_cell->atoms;
		 delete crystal_cell;
		 delete[] crystal_cuts->planes;
		 delete crystal_cuts;
		return 0;
	};

	// cut crystal_cell
	crystal_cut = *crystal_cell;
	for (int i=0; i<crystal_cuts->number; i++)
	{
		// already cut crystal becomes crystal to cut
		crystal_to_cut = crystal_cut;
		// cut it
		crystal_cut = slash(crystal_to_cut, crystal_cuts->planes[i]);
		// delete previous crystal
		delete[] crystal_to_cut.atoms;
		// first iteration deletes crystal_cell->atoms
	};

	if (crystal_cut.number > 0)
	{	// atoms left after cuts
		tofile(argv[2], &crystal_cut); //writes data to file
	} else
		cout << "\nNo atoms left after cutting crystal.\nNo file written\n";

	//free memory
	 delete[] crystal_cut.atoms;
	 delete crystal_cell;
	 delete[] crystal_cuts->planes;
	 delete crystal_cuts;
	 
	return 0;
}

//---------------------------------------------------------------------------

