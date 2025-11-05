#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>

using namespace std;

void max(int n, double *v, int *imax, double *vmax)
{
  *imax = 0;
  *vmax = v[0];
  for(int i= 1; i<n; i++)
    if(v[i] > *vmax)
    {
      *vmax = v[i];
      *imax = i;
    }
}

void min(int n, double *v, int *imin, double *vmin)
{
  *imin = 0;
  *vmin = v[0];
  for(int i= 1; i<n; i++)
    if(v[i] < *vmin)
    {
      *vmin = v[i];
      *imin = i;
    }
}

double mat_avg(int n, int ix, double **mat)
{
  // verificar se eh ix,j ou j,ix
  double avg = 0.0;
  for(int j=0; j<n; j++)
    avg += mat[j][ix];
  return avg/(double)n;
}

int main(int argc, char **argv)
{

  if (argc < 2)
  {
    std::cout << "\n Usage: calc_def <base_dir>\n" << std::endl;
    return 1;
  }

  int n = 4000, nx = 61, ny = 61;
  double value;
  double v[4000];
  std::string base = argv[1];

  double **matc = new double* [nx];
  for(int i=0; i<nx; i++)
    matc[i] = new double[ny];

  for(int i=0; i<nx; i++)
    for(int j=0; j<ny; j++)
      matc[i][j] = 0.0;

  // array para guardar as componentes de F
  static double strain[4000][61][61];

  for(int k=0; k<n; k++)
  {
    std::stringstream stream;
    stream << base;
    stream << "/fd_output_000";
    stream << std::setfill('0') << std::setw(4) << k;
    stream << ".dat";
    std::string fname = stream.str();

    cout << "Lendo " << fname << endl;

    // Le todos arquivos
    int numnodes;
    ifstream ifile;
    ifile.open(fname.c_str());
    ifile >> numnodes;
    for(int i=0; i<nx; i++)
    {
      for(int j=0; j<ny; j++)
      {
        double f11,f12,f13,f21,f22,f23,f31,f32,f33;
        ifile >> f11 >> f12 >> f13;
        ifile >> f21 >> f22 >> f23;
        ifile >> f31 >> f32 >> f33;

	// Calcula C11
	strain[k][i][j] = 0.5 * ((f11*f11 + f21*f21) - 1.0);
	
	// Calcula C22
	//strain[k][i][j] = 0.5 * ((f12*f12 + f22*f22) - 1.0);	

      }
    }
    ifile.close();
  }

  // Calcula max de C11 por nó
  double maxc11;
  for(int i=0; i<nx; i++)
  {
    for(int j=0; j<ny; j++)
    {
      maxc11 = strain[0][i][j];
      for(int k=1; k<n; k++)
      {
	if (strain[k][i][j] > maxc11)
	  maxc11 = strain[k][i][j];
      }
      matc[i][j] = maxc11;
    }
  }
  
  // Escreve max C11 no tempo para cada noh em arquivo
  ofstream ofile;
  ofile.open("maxdef.dat");
  for(int i=0; i<nx; i++)
  {
    for(int j=0; j<ny; j++)
    {
      ofile << scientific << matc[i][j] << " ";
    }
    ofile << endl;
  }
  ofile.close();

  // Escreve arquivo para plotar X vs avg def C11 na direcao y

  double **temp = new double* [nx];
  for(int i=0; i<nx; i++)
    temp[i] = new double[ny];

  for(int i=0; i<61; i++)
    for(int j=0; j<61; j++)
      temp[i][j] = strain[1000][i][j];

  ofstream ofileavg;
  ofileavg.open("avgdef.dat");
  for(int i=0; i<nx; i++)
  {
    double avgy = mat_avg(ny,i,matc);
    //double avgy = mat_avg(ny,i,temp);
    double xdst = (double)i/nx;
    ofileavg << scientific << xdst << " " << avgy;
    ofileavg << endl;
  }
  ofileavg.close();

  // Free memory
  for(int i=0; i<nx; i++)
    delete [] matc[i];
  delete [] matc;

  return 0;
}
