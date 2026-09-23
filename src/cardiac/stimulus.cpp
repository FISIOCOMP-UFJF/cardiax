#include <stdexcept>
#include "stimulus.hpp"
#include "util/pugixml.hpp"

void Stimuli::add(double start, double duration, double sval,
                  double x0, double x1, double y0, double y1,
                  double z0, double z1)
{
  Stimulus *s = new Stimulus(start, duration, sval, x0, x1, y0, y1, z0, z1);
  svec.push_back(s);
  num_stimuli++;
}

void Stimuli::clear()
{
  for(uint i=0; i<svec.size(); i++ ) {
    delete svec[i];
  }
  svec.clear();
}

// Check and apply for Monodomain model
void Stimuli::check(const double time, const Mesh &msh,
                    std::set<uint> & snodes, double * sval, bool * apply)
{
  bool stim, sregion;
  arma::vec3 pt;
  std::vector<Stimulus*>::iterator p;
  const std::vector<arma::vec3> & pts = msh.get_points();

  *apply = false;

  if(num_stimuli != 0)
  {    
    for(p = svec.begin()+num_applied; p != svec.end(); ++p)
    {
      stim = (time >= (*p)->start) && (time <= (*p)->start + (*p)->duration);
      if (stim)
      {
        Stimulus *s = *p;
        *sval = s->value;
        *apply = true;

        //cout << " Applying stimulus at time " << time << endl;
        for(uint i=0; i<pts.size(); i++)
        {
          pt = pts[i];

          // stimulate a cuboid region
          sregion  = (pt(0) >= s->x0) && (pt(0) <= s->x1);
          sregion &= (pt(1) >= s->y0) && (pt(1) <= s->y1);
          sregion &= (pt(2) >= s->z0) && (pt(2) <= s->z1);
          
          // stimulate an ellipsoidal region                  
          //double X1,X2,X3;
          //double A1=0, A2=0,A3=0;
          //double R1=0.5*20000, R2=0.5*20000, R3=0.5*43000;
          //double A1=0,    A2=0,    A3=-6900;
          //double R1=8000, R2=8000, R3=11500;
          //X1=pt(0); X2=pt(1); X3=pt(2);
          //sregion = ( (X1-A1)*(X1-A1)/(R1*R1) +
          //            (X2-A2)*(X2-A2)/(R2*R2) +
          //            (X3-A3)*(X3-A3)/(R3*R3) ) <= 1.0;

          if (sregion) snodes.insert(i);
        }
      }

      if ( time >= ((*p)->start + (*p)->duration) )
      {
        num_applied++;
        num_stimuli--;

        // remove the recently applied stimulus
        //svec.erase(p);
      }

    }
  }

}

//
// Bidomain model only
//
void Stimuli::check_bidomain(const double time, const Mesh &msh,
                             std::map<uint,uint> & tmap,
                             std::set<uint> & snodes,
                             double * sval, bool * apply)
{
  bool stim, sregion;
  arma::vec3 pt;
  std::map<uint,uint>::iterator it;
  std::vector<Stimulus*>::iterator p;

  const std::vector<arma::vec3> & pts = msh.get_points();

  *apply = false;

  if(num_stimuli != 0)
  {
    for(p = svec.begin()+num_applied; p != svec.end(); ++p)
    {
      stim = (time >= (*p)->start) && (time <= (*p)->start + (*p)->duration);
      if (stim)
      {
        Stimulus *s = *p;
        *sval = s->value;
        *apply = true;

        for(uint i=0; i<pts.size(); i++)
        {
          pt = pts[i];

          // cuboid region for stimulus
          sregion  = (pt(0) >= s->x0) && (pt(0) <= s->x1);
          sregion &= (pt(1) >= s->y0) && (pt(1) <= s->y1);
          sregion &= (pt(2) >= s->z0) && (pt(2) <= s->z1);

          // stimulate an ellipsoidal region
          //a1, a2 = 28000, 43000
          //c1, c2 = 64000, 70000
          //double X1,X2,X3,A1=0,A2=0,A3=-15000,R1=34000,R2=34000,R3=50000;
          //X1=pt(0); X2=pt(1); X3=pt(2);
          //sregion = ( (X1-A1)*(X1-A1)/(R1*R1) +
          //            (X2-A2)*(X2-A2)/(R2*R2) +
          //            (X3-A3)*(X3-A3)/(R3*R3) ) <= 1.0;

          //double X1,X2,X3;
          //double A1=0    , A2=0,   A3=0;
          //double R1=20000, R2=20000, R3=44000;
          //X1=pt(0); X2=pt(1); X3=pt(2);
          //sregion = ( (X1-A1)*(X1-A1)/(R1*R1) +
          //            (X2-A2)*(X2-A2)/(R2*R2) +
          //            (X3-A3)*(X3-A3)/(R3*R3) ) <= 1.0;

          if (sregion)
          {
            it = tmap.find(i);
            if (it != tmap.end() )
              snodes.insert(it->second);
          }
        }

        if ( time >= ((*p)->start + (*p)->duration) )
        {
          num_applied++;
          num_stimuli--;
        }

      }
    }
  }
}

void Stimuli::read(const string & filename)
{

  double start, duration;
  double x0, x1, y0, y1, z0, z1, strength;
  string sfile = filename + ".stim";

  ifstream in;
  in.open(sfile.c_str());
  in >> num_stimuli;

  svec.reserve(num_stimuli);

  cout << "Number of stimuli: " << num_stimuli << endl;

  for(int i=0; i<num_stimuli; i++)
  {
    in >> start >> duration >> strength >> x0 >> x1 >> y0 >> y1 >> z0 >> z1;
    Stimulus *s = new Stimulus(start, duration, strength, x0, x1, y0, y1, z0, z1);
    svec.push_back(s);
    cout << "Stimulus " << i+1 << " : ";
    cout << "[" << start << "," << start + duration << "]" << endl;
  }

  in.close();
}

void Stimuli::read_xml(const string & filename)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(filename.c_str());

  if (!result)
  {
    std::cout << "XML parsed with errors, attr value: ["
              << doc.child("mesh").attribute("attr").value() << "]\n";
    std::cout << "Error description: " << result.description() << "\n";
    std::cout << "Error offset: " << result.offset;
    std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
  }

  // element type
  pugi::xml_node stims = doc.child("electrophysiology").child("stimuli");

  int i = 0;
  num_stimuli = stims.attribute("number").as_int();
  svec.reserve(num_stimuli);
  cout << "Number of stimuli: " << num_stimuli << endl;
  for (pugi::xml_node_iterator it = stims.begin(); it != stims.end(); ++it)
  {
    double start = it->attribute("start").as_double();
    double duration = it->attribute("duration").as_double();
    double value = it->attribute("value").as_double();
    double x0 = it->attribute("x0").as_double();
    double x1 = it->attribute("x1").as_double();
    double y0 = it->attribute("y0").as_double();
    double y1 = it->attribute("y1").as_double();
    double z0 = it->attribute("z0").as_double();
    double z1 = it->attribute("z1").as_double();
    Stimulus *s = new Stimulus(start, duration, value, x0, x1, y0, y1, z0, z1);
    svec.push_back(s);
    cout << " Stimulus " << i+1 << ": [" << start
         << "," << start + duration
         << ", value: " << value << " ]" << endl;
    i++;
  }

}


int Stimuli::read_toml(const toml::table & cfg, const std::string & key)
{
  const toml::array * regions = cfg.at_path(key + ".regions").as_array();

  if (!regions || regions->empty())
  {
    cout << "Number of stimuli (TOML " << key << "): 0" << endl;
    return 0;
  }

  svec.reserve(svec.size() + regions->size());

  int i = 0;
  for (const toml::node & node : *regions)
  {
    i++;
    const std::string where = key + ".regions entry " + std::to_string(i);

    const toml::table * t = node.as_table();
    if (!t)
      throw std::runtime_error(where + ": must be an inline table { ... }");

    auto scalar = [&](const char * name) -> double {
      auto v = (*t)[name].value<double>();          // also accepts integers
      if (!v) throw std::runtime_error(where + ": missing or invalid '" + name + "'");
      return *v;
    };

    auto vec3 = [&](const char * name, double out[3]) {
      const toml::array * a = (*t)[name].as_array();
      if (!a || a->size() != 3)
        throw std::runtime_error(where + ": '" + name + "' must be [x, y, z]");
      for (int d = 0; d < 3; d++)
      {
        auto v = (*a)[d].value<double>();
        if (!v) throw std::runtime_error(where + ": non-numeric value in '" + name + "'");
        out[d] = *v;
      }
    };

    const double start    = scalar("start");
    const double duration = scalar("duration");
    const double value    = scalar("value");
    double lo[3], hi[3];
    vec3("min", lo);
    vec3("max", hi);

    if (duration <= 0.0)
      throw std::runtime_error(where + ": 'duration' must be > 0");
    for (int d = 0; d < 3; d++)
      if (lo[d] > hi[d])
        throw std::runtime_error(where + ": min > max in component " + std::to_string(d));

    svec.push_back(new Stimulus(start, duration, value,
                                lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]));

    cout << " Stimulus " << i << ": [" << start << "," << start + duration
         << "] value " << value
         << "  box [" << lo[0] << "," << hi[0] << "]x[" << lo[1] << "," << hi[1]
         << "]x[" << lo[2] << "," << hi[2] << "]" << endl;
  }

  num_stimuli = (int) svec.size();
  cout << "Number of stimuli (TOML " << key << "): " << i << endl;
  return i;
}