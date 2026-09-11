#include "eikonal_reactive.hpp"
#include "mesh/writer_hdf5.hpp"
#include <iostream>
#include <cmath>
#include "util/pugixml.hpp"


using namespace std;

EikonalReactive::EikonalReactive() : CardiacProblem(), stim_apply_nodes(false) {
    mesh = new Mesh();
    writer = new WriterHDF5(mesh);

    parameters.rename("Eikonal_parameters");
    parameters.add("vel_f", 0.006);
    parameters.add("vel_s", 0.0002);
    parameters.add("vel_n", 0.0002);
}

EikonalReactive::~EikonalReactive() {
    delete cells;
    delete cellmodel;
}

void EikonalReactive::setup(std::string & b, std::string & c, std::string & m,
                            double dt, double T, double pr, double pa) {
    timestep  = dt;
    totaltime = T;
    printrate = pr;
    printrate_apd = pa;
  
    mesh_filename = b;
    stimuli_filename = b;

    cell_name = c;
    odesolver = m;
}

void EikonalReactive::init() {
    tip = TimeParameters(timestep, totaltime, printrate);

    mesh->read_xml(mesh_filename);
    stimuli.read_xml(stimuli_filename);

    fespace.set_mesh(mesh);
    ndofs = mesh->get_n_points();

    stim_values.zeros(ndofs);

    std::size_t pos  = mesh_filename.find(".xml");
    std::string output = mesh_filename.substr(0,pos) + "_output";
    int nsteps = tip.get_size(); 
    writer->open(output, nsteps+1, timestep);

    eikonal_solver.set_velocities(parameters["vel_f"], parameters["vel_s"], parameters["vel_n"]);
    
    eikonal_solver.solve(mesh, mesh_filename);
    
    
    lat = eikonal_solver.get_lat();

    cellmodel = CellModel::create(cell_name);
    cellmodel->setup(odesolver, timestep, totaltime, 1.0);
    cells = new Cells(ndofs, cellmodel);
    cells->init();
}

void EikonalReactive::read_eikonal_solution(const std::string &mshfile) {
    std::cout << " -- Initializing local activation time (Reading from file) --" << std::endl; 
    
    pugi::xml_document doc;
    doc.load_file(mshfile.c_str());

    lat.zeros(ndofs);
    pugi::xml_node eikonal_data = doc.child("mesh").child("eikonal");
    int n_read = 0;
    
    if (eikonal_data.child("node")) {
        for(pugi::xml_node node = eikonal_data.child("node"); node; node = node.next_sibling("node")) {
            int index = node.attribute("id").as_int(); 
            if(index >= 0 && index < (int)ndofs) {
                // The per-node LAT is given in MILLISECONDS.
                lat(index) = node.attribute("lat").as_double(); 
                n_read++;
            }
        }
    }
    
    if(n_read != (int) ndofs) {
        pugi::xml_node pvloop_data = doc.child("mesh").child("pvloop");
        double begin_active_stress = 0.0;
        
        if(pvloop_data) {
            begin_active_stress = pvloop_data.attribute("passive_time").as_double() * ms_to_solver_time(); 
        }
            
        lat.fill(begin_active_stress);
        std::cout << " No valid/complete per-node LAT found. Using uniform activation: " 
                  << begin_active_stress << std::endl;
    }
}

void EikonalReactive::initial_conditions() {
    tip.reset();
    cells->init();

    warmup_cells(lat);

    int lat_idx = cells->get_model().lat_var_index();
    if (lat_idx >= 0)
        cells->set_var(lat_idx, lat);

    cells->advance(tip.time(), timestep, stim_val, stim_nodes);
}

void EikonalReactive::advance() {
    if( !tip.finished() ) {
        tip.increase_time();
        timer.enter("ODEs");
        solve_odes();
        timer.leave();
    }
}

void EikonalReactive::solve_odes() {
    stimuli.check(tip.time(), *mesh, stim_nodes, &stim_val, &stim_apply);
  
    if (stim_apply) {
        cells->advance(tip.time(), timestep, stim_val, stim_nodes);
        stim_nodes.clear();
    } else if(stim_apply_nodes) {
        cells->advance(tip.time(), timestep, stim_values);
        stim_values.fill(0);
        stim_apply_nodes = false;
    } else {
        cells->advance(tip.time(), timestep);
    }
}

void EikonalReactive::set_stimulus_value(int index, double val) {
    stim_apply_nodes = true;
    stim_values(index) = val;
}

void EikonalReactive::apply_lat_stimulus(double amplitude, double duration, double period) {
    if (stim_values.n_elem != lat.n_elem) {
        stim_values.zeros(lat.n_elem);
    }

    double t = tip.time();
    if (period > 0.0) t = std::fmod(t, period);

    bool any = false;
    for (arma::uword i = 0; i < lat.n_elem; i++) {
        if (t >= lat(i) && t < lat(i) + duration) {
            stim_values(i) = amplitude;
            any = true;
        }
    }

    if (any) stim_apply_nodes = true;
}

void EikonalReactive::solve() {
    //Do nothing
}

void EikonalReactive::set_conductivity(int cond) {
    eikonal_solver.set_conductivity(cond);
}