#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>
#include <stdlib.h>
#include "apd.h"

using namespace std;

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w=50)
{
    if ( (x != n) && (x % (n/100) != 0) ) return;
 
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

int main(int argc, char **argv)
{
  if (argc < 5)
  {
    cout << "\n Usage: calc_apd <base_dir> <n_files> <n_nodes> <dt>\n" << endl;
    return 1;
  }

  int num_files, num_nodes;
  double timestep, value;

  string base = argv[1];
  num_files = atoi(argv[2]);
  num_nodes = atoi(argv[3]);
  timestep = atof(argv[4]);

  double *mapd, *actv, *repo, *v, **vm;
  mapd = new double[num_nodes];
  actv = new double[num_nodes];
  repo = new double[num_nodes];
  v    = new double[num_files];  
  vm   = new double*[num_files];
  for(int i=0; i<num_files; i++)
    vm[i] = new double[num_nodes];

  cout << "Reading files" << endl;  
  for(int k=0; k<num_files; k++)
  {
    stringstream stream;
    stream << base;
    stream << "vm_";
    stream << setfill('0') << setw(5) << k;
    stream << ".txt";
    string fname = stream.str();

    loadbar(k, num_files);

    ifstream ifile;
    ifile.open(fname.c_str());
    if(!ifile)
    {
      cerr << "Erro lendo arquivo " << fname << endl;
      return 1;
    }
    for(int i=0; i<num_nodes; i++)
    {
      ifile >> value;
      vm[k][i] = value;
    }
    ifile.close();
  }

  cout << "Calculating APDfor each node" << endl;
  for(int i=0; i<num_nodes; i++)
  {
    for(int k=0; k<num_files; k++) v[k] = vm[k][i];
    double lapd, lact, lrep;
    
    calc_apd(timestep, num_files, v, 90.0, &lapd, &lact, &lrep);
    mapd[i] = lapd;
    actv[i] = lact;
    repo[i] = lrep;
    
    cout << "Node " << i << " APD90 " << mapd[i] << endl;
  }

  cout << "Saving files" << endl;
  string outfile1 = base + "mapcolorAPD.dat";
  ofstream ofile1(outfile1.c_str());
  for(int i=0; i<num_nodes; i++)
    ofile1 << scientific << mapd[i] << endl;
  ofile1.close();  
  cout << "APD data written in file: " << outfile1 << endl;
  
  string outfile2 = base + "mapcolorACT.dat";
  ofstream ofile2(outfile2.c_str());
  for(int i=0; i<num_nodes; i++)
    ofile2 << scientific << actv[i] << endl;
  ofile2.close();  
  cout << "ACT data written in file: " << outfile2 << endl;
  
  string outfile3 = base + "mapcolorREP.dat";
  ofstream ofile3(outfile3.c_str());
  for(int i=0; i<num_nodes; i++)
    ofile3 << scientific << repo[i] << endl;
  ofile3.close();  
  cout << "REP data written in file: " << outfile3 << endl;
  
  // free memory
  delete [] v;
  delete [] mapd;
  delete [] actv;
  delete [] repo;
  for(int i=0; i<num_files; i++)
    delete [] vm[i];
  delete [] vm;

  cout << "Done" << endl;

  return 0;
}

