#include "cardiacproblem.hpp"

CardiacProblem::CardiacProblem()
{
  // setup default parameter values
  parameters.add("surface_to_volume", 0.14);
  parameters.add("theta_method", 0.5);
  parameters.add("pcgtol", 1.0e-16);
  parameters.add("maxnz", 40);
}

CardiacProblem::~CardiacProblem()
{
  writer->close();
  delete mesh;
  delete writer;
}

void CardiacProblem::setup(std::string & b, std::string & c, std::string & m,
                           double dt, double T, double pr, double pa)
{
  cout << "Setup of cardiac problem" << endl;
  const double theta = parameters["theta_method"];

  timestep  = dt;
  totaltime = T;
  printrate = pr;
  printrate_apd = pa;
  
  mesh_filename = b;
  stimuli_filename = b;

  cell_name = c;
  odesolver = m;

  cout << "Parabolic solver: ";
  if (theta == 0.0)      cout << "Explicit Euler" << endl;
  else if (theta == 1.0) cout << "Implicit Euler" << endl;
  else if (theta == 0.5) cout << "Crank-Nicolson" << endl;
}

void CardiacProblem::setup_types(std::string & f)
{
  if( !file_exists(f) )
    cout << "Cells: all cells are of the same type" << endl;
  else
  {
    cout << "Cells: configuring cell types" << endl;
   
    int aux, size;
    int * vtypes;
    ifstream in(f.c_str());
    in >> size;
    vtypes = new int[size];

    for(int i=0; i<size; i++)
    {
      in >> aux;
      vtypes[i] = aux;
    }

    cells->set_cell_types(size, vtypes);   
    delete [] vtypes;
  }
}

void CardiacProblem::write_data(const arma::vec & u, const std::string & s, int * step)
{
  if (tip.time2print())
  {
    std::string aux = s.substr(s.length() - 2);

    if (aux == "ve")
      writer->write_ve_step(*step, u.memptr());
    else if (aux == "vm")
      writer->write_vm_step(*step, u.memptr());
    else if(aux=="dv")
    {
      arma::vec current = u( arma::span(u.size()/2, u.size()-1) );
      writer->write_ve_step(*step, current.memptr());
    }
    
    *step = *step + 1;
  }
}

void CardiacProblem::write_data(const arma::vec & u, const arma::vec & displ,
                                const std::string & s, int step)
{
    if (tip.time2print())
    {
        std::string aux = s.substr(s.length()-2);

        // potentials scalar field
        if(aux=="ve")
          writer->write_ve_step(step, u.memptr());
        else if(aux=="vm")
            writer->write_vm_step(step, u.memptr());

        // displacement vector field
        writer->write_displ_step(step, displ.memptr());

  }
}

void CardiacProblem::write_data_text(const arma::vec & vm, int * step)
{
  //if (tip.time2print())
  //{
    char buffer[256];
    sprintf(buffer,"output/vm_text/vm_%05d.txt", *step);
    string name(buffer);   
    vm.save(name, arma::raw_ascii);
    *step = *step + 1;
  //}
}
