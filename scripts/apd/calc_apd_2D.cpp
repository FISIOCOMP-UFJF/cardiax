#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>
#include "apd.h"

using namespace std;

/**
   Simple program to compute the APD of the files generated from
   my old electromechanical simulator (bidopetsc3 + Flagshyp).
*/

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    std::cout << "\n Usage: calc_apd_2D <base_dir>\n" << std::endl;
    return 1;
  }

  int n = 4000, nx = 61, ny = 61;
  double value;
  double v[4000];
  std::string base = argv[1];

  double **mapd = new double* [61];
  for(int i=0; i<61; i++)
    mapd[i] = new double[61];

  static double vm[4000][61][61];

  for(int k=0; k<n; k++)
  {
    std::stringstream stream;
    stream << base;
    stream << "/vm_text_000";
    stream << std::setfill('0') << std::setw(4) << k;
    stream << ".txt";
    std::string fname = stream.str();

    cout << "Lendo " << fname << endl;

    // Le todos arquivos
    ifstream ifile;
    ifile.open(fname.c_str());
    for(int i=0; i<nx; i++)
    {
      for(int j=0; j<ny; j++)
      {
	ifile >> value;
	vm[k][i][j] = value;
	//Vm(k,i,j) = value;
      }
    }
    ifile.close();
  }

  // Calcula APD por n�
  for(int i=0; i<nx; i++)
  {
    for(int j=0; j<ny; j++)
    {
      for(int k=0; k<n; k++) v[k] = vm[k][i][j];

      double lapd, lact, lrep;
      calc_apd(0.1,n,v,90.0, &lapd, &lact, &lrep);
      mapd[i][j] = lapd;
      
      //mapd[i][j] = calc_apd(0.1,n,v,90.0);
      //mapd[i][j] = calc_rep_time(0.1,n,v,90.0);

      cout << " Calculando APD90 noh (";
      cout << i << "," << j << ") ";
      cout << mapd[i][j] << endl;
    }
  }

  ofstream ofile;
  ofile.open("mapcolor_apd.dat");
  for(int i=0; i<nx; i++)
  {
    for(int j=0; j<ny; j++)
    {
      ofile << scientific << mapd[i][j] << " ";
    }
    ofile << endl;
  }
  ofile.close();

  // Plot x vs APD_avg in y direction
  ofstream ofileavg;
  ofileavg.open("apdavg_y_vs_x.dat");
  for(int i=0; i<nx; i++)
  {
    double avgy = apd_avg(ny,i,mapd);
    double xdst = (double)i/nx;
    ofileavg << scientific << xdst << " " << avgy;
    ofileavg << endl;
  }
  ofileavg.close();

  // Free memory
  for(int i=0; i<nx; i++)
    delete [] mapd[i];
  delete [] mapd;

  return 0;
}
