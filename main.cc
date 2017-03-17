#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
int number_of_unit_cells;
int number_of_atoms;
int n_cell = 0;
int n;
double potential;
double constant = 5.26e-10; //reduced unit 5.26 angstrom
double eps = 0.0104;
double sigma = 3.4e-10;
double atom_position[100000][4][3];
double atoms[100000][3];
double a = 5.26e-10;

//Create the FCC Lattice
int init (int number_of_unit_cells)
{
  int x, y, z, i;
  i = 0;
  double half = constant / 2;

  for (x = 0; x < number_of_unit_cells; x++)
  {
    double xstep = constant * x;
    for (y = 0; y < number_of_unit_cells; y++)
    {
      double ystep = constant * y;
      for (z = 0; z < number_of_unit_cells; z++)
      {
        double zstep = constant * z;

  
    atom_position[i][0][0] = xstep;
    atom_position[i][0][1] = ystep;
    atom_position[i][0][2] = zstep;

    atom_position[i][1][0] = xstep + half;
    atom_position[i][1][1] = ystep + half;
    atom_position[i][1][2] = zstep;

    atom_position[i][2][0] = xstep + half;
    atom_position[i][2][1] = ystep;
    atom_position[i][2][2] = zstep + half;

    atom_position[i][3][0] = xstep;
    atom_position[i][3][1] = ystep + half;
    atom_position[i][3][2] = zstep + half;
    i++;  
      }
    }
  }
  return 0;
}

//Calculate LJ potential energy with no periodic boundary conditions
double energy (int p)
{
int counter = 0;
double distance;
potential = 0;
for (int i = 0; i < p - 1; i++)
{
  for (int j = i + 1; j < p; j++)
  {
  distance = sqrt(pow(atoms[i][0] - atoms[j][0], 2) + pow(atoms[i][1] - atoms[j][1], 2) + pow(atoms[i][2] - atoms[j][2], 2)) * 1 / sigma;
    potential = potential + 4 * (pow(distance, -12.) - pow(distance, -6.));
    counter++;
  }
}
return potential / p;
}

//Calculate LJ potential energy with periodic boundary conditions
double energy_pbc (int n, int p)
{
int counter = 0;
double distance;
double xij, yij, zij;
double epotential = 0;
double fraction = 0.5;
for (int i = 0; i < p - 1; i++)
{
  for (int j = i + 1; j < p; j++)
  {
  xij = atoms[i][0] - atoms[j][0];
  yij = atoms[i][1] - atoms[j][1];
  zij = atoms[i][2] - atoms[j][2];
  xij = xij / a;
  yij = yij / a;
  zij = zij / a;

  if (xij > n * fraction)
    xij = xij - n;
  else if (xij < - n * fraction)
    xij = xij + n;

  if (yij > n * fraction)
    yij = yij - n;
  else if (yij < - n * fraction)
    yij = yij + n;

  if (zij > n * fraction)
    zij = zij - n;
  else if (zij < - n * fraction)
    zij = zij + n;


    distance = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2)) * a / sigma;
    epotential = epotential + 4 * (pow(distance, -12.) - pow(distance, -6.));
    counter++;
  }
}
return epotential / p;
}

 
//Write the XYZ with the initial coordinates of the system and create a new array easier to work with. Return the total number of atoms.
int write_text ()
{
  ofstream myfile;
  myfile.open("coordinates.xyz");
  myfile << 4 * n * n * n << endl;
  myfile << endl;
  int p=0;
  for (int i = 0; i < n * n * n; i++){
    atoms[p][0] = atom_position[i][0][0];        
    atoms[p][1] = atom_position[i][0][1];        
    atoms[p][2] = atom_position[i][0][2];        
    p++;
    atoms[p][0] = atom_position[i][1][0];        
    atoms[p][1] = atom_position[i][1][1];        
    atoms[p][2] = atom_position[i][1][2];        
    p++;
    atoms[p][0] = atom_position[i][2][0];        
    atoms[p][1] = atom_position[i][2][1];        
    atoms[p][2] = atom_position[i][2][2];        
    p++;
    atoms[p][0] = atom_position[i][3][0];        
    atoms[p][1] = atom_position[i][3][1];        
    atoms[p][2] = atom_position[i][3][2];        
    p++;
    myfile << "Al " << atom_position[i][0][0] << " " << atom_position[i][0][1] << " " << atom_position[i][0][2] << endl;
    myfile << "Al " << atom_position[i][1][0] << " " << atom_position[i][1][1] << " " << atom_position[i][1][2] << endl;
    myfile << "Al " << atom_position[i][2][0] << " " << atom_position[i][2][1] << " " << atom_position[i][2][2] << endl;
    myfile << "Al " << atom_position[i][3][0] << " " << atom_position[i][3][1] << " " << atom_position[i][3][2] << endl;
  }
  return p;
}

  

int main()
{
  cout << "give number of unit cells" << endl;
  cin >> n;
//  cout << "lattice constant" << endl;
//  cin >> constant;
  init(n);
  number_of_atoms = write_text();
  printf("%d \n", number_of_atoms);
  cout << energy(number_of_atoms) * eps << endl;
  cout << energy_pbc(n, number_of_atoms) * eps << endl;
}
