#include "mesh.hpp"
#include "fem/fe.hpp"
#include "util/pugixml.hpp"

#include <algorithm>

void Mesh::update_point(const int i, const arma::vec3 &pt)
{
  points[i](0) = pt(0);
  points[i](1) = pt(1);
  points[i](2) = pt(2);
}

void Mesh::get_element_pt_nums(int elnum, std::vector<int> &pnums) const
{
  pnums = elems[elnum].get_pt_nums();
}

void Mesh::get_boundary_element_pt_nums(int selnum, std::vector<int> &pnums) const
{
  pnums = boundary_elems[selnum].get_pt_nums();
}

void Mesh::read(const std::string &filename)
{
  cout << "Reading mesh file: " << filename << endl;

  int fibcode;
  string ptsfile, elefile, axefile;
  ifstream inp, ine, inf;

  ptsfile = filename + ".pts";
  elefile = filename + ".elem";
  axefile = filename + ".fib";

  //
  // Check files
  //
  if (!file_exists(ptsfile)) error(".pts file not found");
  if (!file_exists(elefile)) error(".elem file not found");
  if (!file_exists(axefile)) error(".fib file not found");

  //
  // Read points
  //
  inp.open(ptsfile.c_str());
  inp >> n_points;
  for (int i = 0; i < n_points; i++) {
    arma::vec3 p;
    inp >> p[0] >> p[1] >> p[2];
    points.push_back(p);
  }
  inp.close();

  //
  // Elements with fibers
  //
  ine.open(elefile.c_str());
  inf.open(axefile.c_str());
  ine >> n_elements;
  inf >> fibcode;

  cout << "Fiber model: ";
  if (fibcode == 0)
  {
    cout << "orthotropic" << endl;
    proptype = ORTHOTROPIC;
  }
  else
  {
    cout << "transversely isotropic" << endl;
    proptype = TRANSVERSELY_ISOTROPIC;
  }

  for (int i = 0; i < n_elements; i++) {
    int index;
    string elname;
    arma::vec3 fib, trn, nor;

    ine >> elname;

    if (elname == "Ln") {
      n_dim = 1;
      n_node = 2;
      std::vector<int> pnum(2);
      ine >> pnum[0] >> pnum[1] >> index;
      read_fiber(inf, fibcode, fib, trn, nor);
      elems.push_back(Element(ELEM_SEGM, pnum, fib, trn, nor, index));
    } else if (elname == "Qd") {
      n_dim = 2;
      n_node = 4;
      std::vector<int> pnum(4);
      ine >> pnum[0] >> pnum[1] >> pnum[2] >> pnum[3] >> index;
      read_fiber(inf, fibcode, fib, trn, nor);
      elems.push_back(Element(ELEM_QUAD, pnum, fib, trn, nor, index));
    } else if (elname == "Tr") {
      n_dim = 2;
      n_node = 3;
      std::vector<int> pnum(3);
      ine >> pnum[0] >> pnum[1] >> pnum[2] >> index;
      read_fiber(inf, fibcode, fib, trn, nor);
      elems.push_back(Element(ELEM_TRIG, pnum, fib, trn, nor, index));
    } else if (elname == "Tt") {
      n_dim = 3;
      n_node = 4;
      std::vector<int> pnum(4);
      ine >> pnum[0] >> pnum[1] >> pnum[2] >> pnum[3] >> index;
      read_fiber(inf, fibcode, fib, trn, nor);
      elems.push_back(Element(ELEM_TETRA, pnum, fib, trn, nor, index));
    } else if (elname == "Hx") {
      n_dim = 3;
      n_node = 8;
      std::vector<int> pnum(8);
      ine >> pnum[0] >> pnum[1] >> pnum[2] >> pnum[3]
          >> pnum[4] >> pnum[5] >> pnum[6] >> pnum[7] >> index;
      read_fiber(inf, fibcode, fib, trn, nor);
      elems.push_back(Element(ELEM_HEXA, pnum, fib, trn, nor, index));
    } else
      print_error("mesh.cpp", "reading", "element not found");
  }
  ine.close();
  inf.close();

  //
  // Surface elements
  //
  string sfile = filename + ".surf";
  ifstream ise(sfile.c_str());
  if (ise) {
    // Read surface elements
    ise >> n_boundary_elements;
    boundary_elems.reserve(n_boundary_elements);
    string elname;
    for (int i = 0; i < n_boundary_elements; i++) {
      ise >> elname;
      if (elname == "Ln") {
        int index;
        std::vector<int> pnum(2);
        ise >> pnum[0] >> pnum[1] >> index;
        boundary_elems.push_back(Element(ELEM_SEGM, pnum, index));
      } else
        print_error("mesh.cpp", "reading", "surface element not found in mesh");
    }
    ise.close();
  }

  //
  // Constrained boundary nodes
  //
  string bfile = filename + ".bdn";
  ifstream ibd(bfile.c_str());
  int num_boundary_nodes, node;
  double value;
  if (ibd) {
    // Read constrained nodes
    ibd >> num_boundary_nodes;
    for (int i = 0; i < num_boundary_nodes; i++) {
      ibd >> node >> value;
      bnodes_map[node] = value;
    }
    ibd.close();
  }
}

void Mesh::read_xml(const std::string &filename)
{
  std::string extension = file_extension(filename);
  if (extension == "gz") {
    cout << "GZipped XML file: please uncompress with gunzip." << endl;
    exit(0);
  }

  cout << "Reading XML mesh file: " << filename << endl;
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(filename.c_str());

  if (!result) {
    std::cout << "XML parsed with errors, attr value: ["
              << doc.child("mesh").attribute("attr").value() << "]\n";
    std::cout << "Error description: " << result.description() << "\n";
    std::cout << "Error offset: " << result.offset;
    std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
  }

  // element type
  std::string celltype = doc.child("mesh").attribute("celltype").value();

  // number of dimensions
  n_dim = doc.child("mesh").attribute("dim").as_int();

  //
  // reading nodes
  //
  pugi::xml_node nodes = doc.child("mesh").child("nodes");
  for (pugi::xml_node node = nodes.child("node"); node;
       node = node.next_sibling("node")) {
    arma::vec3 p;
    p(0) = node.attribute("x").as_double()*1e1;
    p(1) = node.attribute("y").as_double()*1e1;
    p(2) = node.attribute("z").as_double()*1e1;
    points.push_back(p);
  }

  //
  // reading elements
  //
  pugi::xml_node elements = doc.child("mesh").child("elements");

  if (celltype == "line") {
    n_node = 2;
    for (pugi::xml_node elem = elements.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(2);
      int marker = elem.attribute("marker") ? elem.attribute("marker").as_int() : 0;
      int aha = elem.attribute("aha") ? elem.attribute("aha").as_int() : 0;
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      elems.push_back(Element(ELEM_SEGM, pnum, marker, aha));
    }
  }

  if (celltype == "triangle") {
    n_node = 3;
    for (pugi::xml_node elem = elements.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(3);
      int marker = elem.attribute("marker") ? elem.attribute("marker").as_int() : 0;
      int aha = elem.attribute("aha") ? elem.attribute("aha").as_int() : 0;
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      pnum[2] = elem.attribute("v2").as_int();
      elems.push_back(Element(ELEM_TRIG, pnum, marker, aha));
    }
  }

  if (celltype == "quadrilateral") {
    n_node = 4;
    for (pugi::xml_node elem = elements.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(4);
      int marker = elem.attribute("marker") ? elem.attribute("marker").as_int() : 0;
      int aha = elem.attribute("aha") ? elem.attribute("aha").as_int() : 0;
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      pnum[2] = elem.attribute("v2").as_int();
      pnum[3] = elem.attribute("v3").as_int();
      elems.push_back(Element(ELEM_QUAD, pnum, marker, aha));
    }
  }

  if (celltype == "tetrahedron") {
    n_node = 4;
    for (pugi::xml_node elem = elements.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(4);
      int marker = elem.attribute("marker") ? elem.attribute("marker").as_int() : 0;
      int aha = elem.attribute("aha") ? elem.attribute("aha").as_int() : 0;
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      pnum[2] = elem.attribute("v2").as_int();
      pnum[3] = elem.attribute("v3").as_int();
      elems.push_back(Element(ELEM_TETRA, pnum, marker, aha));
    }
  }

  if (celltype == "hexahedron") {
    n_node = 8;
    for (pugi::xml_node elem = elements.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(8);
      int marker = elem.attribute("marker") ? elem.attribute("marker").as_int() : 0;
      int aha = elem.attribute("aha") ? elem.attribute("aha").as_int() : 0;
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      pnum[2] = elem.attribute("v2").as_int();
      pnum[3] = elem.attribute("v3").as_int();
      pnum[4] = elem.attribute("v4").as_int();
      pnum[5] = elem.attribute("v5").as_int();
      pnum[6] = elem.attribute("v6").as_int();
      pnum[7] = elem.attribute("v7").as_int();
      elems.push_back(Element(ELEM_HEXA, pnum, marker, aha));

    //  cout << pnum[0] << " " << pnum[1] << " " << pnum[2] << " " << pnum[3] << " " << pnum[4] << " " << pnum[5] << " " << pnum[6] << " " << pnum[7] << endl;
    }
  }

  if (celltype == "")
    throw std::runtime_error("Error reading XML file: unknown element type");

  // if tri/tet mesh, check for positive volumes
  if (celltype == "triangle") set_tri_positive_volume();
  if (celltype == "tetrahedron") set_tet_positive_volume();

  //
  // reading fibers
  //
  pugi::xml_node fibers = doc.child("mesh").child("element_data");
  if (fibers)
  {
    std::string fiber_model = fibers.attribute("type").value();
    cout << "Fiber model: " << fiber_model << endl;
    if (fiber_model == "fiber_transversely_isotropic")
    {
      proptype = TRANSVERSELY_ISOTROPIC;

      for (pugi::xml_node elem = fibers.child("element"); elem;
           elem = elem.next_sibling("element"))
      {
        int index;
        arma::vec3 f, s, n;
        index = elem.attribute("id").as_int();
        pugi::xml_node fib = elem.child("fiber");
        parse_to_vec3(fib.text().as_string(), f);

        // fill s to be orthogonal to f
        s(0) = -f(1); s(1) = f(0); s(2) = 0.0;
        // compute n
        n = arma::cross(f, s);

        // unit vectors
        f = f / arma::norm(f, 2);
        s = s / arma::norm(s, 2);
        n = n / arma::norm(n, 2);

        if (arma::dot(f,s) > 1e-8 || arma::dot(f,n) > 1e-8 || arma::dot(s,n) > 1e-8)
        {
          cout << "Error!\n";
          cout << "f " << f << endl;
          cout << "s " << s << endl;
          cout << "n " << n << endl;
          cout << "f.s " << arma::dot(f, s) << endl;
          cout << "f.n " << arma::dot(f, n) << endl;
          cout << "s.n " << arma::dot(s, n) << endl;
          exit(0);
        }

        elems[index].set_fiber(f);
        elems[index].set_trans(s);
        elems[index].set_normal(n);
      }
    }
    else if (fiber_model == "fiber_orthotropic")
    {
      proptype = ORTHOTROPIC;

      for (pugi::xml_node elem = fibers.child("element"); elem;
           elem = elem.next_sibling("element"))
      {
        int index;
        arma::vec3 f, s, n;
        index = elem.attribute("id").as_int();
        pugi::xml_node fib = elem.child("fiber");
        pugi::xml_node she = elem.child("sheet");
        pugi::xml_node nor = elem.child("normal");
        parse_to_vec3(fib.text().as_string(), f);
        parse_to_vec3(she.text().as_string(), s);
        parse_to_vec3(nor.text().as_string(), n);
        elems[index].set_fiber(f);
        elems[index].set_trans(s);
        elems[index].set_normal(n);
      }
    }
    //else if (fiber_model == "fiber_isotropic") // isotropic
    else
    {
      proptype = ISOTROPIC;
      for (size_t i = 0; i < elems.size(); i++)
      {
        int index;
        arma::vec3 f, s, n;
        f.zeros();
        s.zeros();
        n.zeros();
        f(0) = 1;
        s(1) = 1;
        n(2) = 1;
        elems[i].set_fiber(f);
        elems[i].set_trans(s);
        elems[i].set_normal(n);
      }
    }
  }


  //
  // reading local basis
  //
  fibers = doc.child("mesh").child("local_basis");
  if (fibers)
  {
    std::string fiber_model = fibers.attribute("type").value();
    cout << "Fiber model: " << fiber_model << endl;
    if (fiber_model == "fiber_orthotropic")
    {

      for (pugi::xml_node elem = fibers.child("element"); elem;
           elem = elem.next_sibling("element"))
      {
        int index;
        arma::vec3 l, c, r;
        index = elem.attribute("id").as_int();
        pugi::xml_node fib = elem.child("long");
        pugi::xml_node she = elem.child("circ");
        pugi::xml_node nor = elem.child("rad");
        parse_to_vec3(fib.text().as_string(), l);
        parse_to_vec3(she.text().as_string(), c);
        parse_to_vec3(nor.text().as_string(), r);
        elems[index].set_long(l);
        elems[index].set_circ(c);
        elems[index].set_rad(r);
      }
    }
    //else if (fiber_model == "fiber_isotropic") // isotropic
    else
    {
      for (size_t i = 0; i < elems.size(); i++)
      {
        int index;
        arma::vec3 f, s, n;
        f.zeros();
        s.zeros();
        n.zeros();
        f(0) = 1;
        s(1) = 1;
        n(2) = 1;
        elems[i].set_long(f);
        elems[i].set_circ(s);
        elems[i].set_rad(n);
      }
    }
  }

  // end reading local basis


  //
  // boundary elements: 2D-edge and 3D-face
  //
  pugi::xml_node bcs = doc.child("mesh").child("boundary");
  std::string bcelltype = bcs.attribute("celltype").value();
  if (bcelltype == "line")
  {
    for (pugi::xml_node elem = bcs.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(2);
      int marker = elem.attribute("marker").as_int();
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      boundary_elems.push_back(Element(ELEM_SEGM, pnum, marker));
    }
  }
  else if (bcelltype == "triangle")
  {
    arma::vec3 ref;
    ref(0) = 0;
    ref(1) = 0;
    ref(2) = 0;
    for (pugi::xml_node elem = bcs.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(3);
      //int idx = elem.attribute("id").as_int();
      int marker = elem.attribute("marker").as_int();
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      pnum[2] = elem.attribute("v2").as_int();
      arma::vec3 A = get_point(pnum[0]);
      arma::vec3 B = get_point(pnum[1]);
      arma::vec3 C = get_point(pnum[2]);

      arma::vec3 AB = (B - A);
      AB = AB / arma::norm(AB, 2);
      arma::vec3 AD = C - A;
      AD = AD / arma::norm(AD, 2);
      arma::vec3 normal = arma::cross(AB, AD);
      arma::vec3 PA = A - ref;
      PA = PA / arma::norm(PA, 2);
/*
      if (arma::dot(PA, normal) < 0)
      {
        int aux = pnum[0];
        pnum[0] = pnum[2];
        pnum[2] = aux;
      }
*/
      boundary_elems.push_back(Element(ELEM_TRIG, pnum, marker));
    }
  }
  else if (bcelltype == "quadrilateral")
  {
    arma::vec3 ref;
    ref(0) = 0;
    ref(1) = 0;
    ref(2) = 0;

    for (pugi::xml_node elem = bcs.child("element"); elem;
         elem = elem.next_sibling("element")) {
      std::vector<int> pnum(4);
      int marker = elem.attribute("marker").as_int();
      pnum[0] = elem.attribute("v0").as_int();
      pnum[1] = elem.attribute("v1").as_int();
      pnum[2] = elem.attribute("v2").as_int();
      pnum[3] = elem.attribute("v3").as_int();

      arma::vec3 A = get_point(pnum[0]);
      arma::vec3 B = get_point(pnum[1]);
      arma::vec3 C = get_point(pnum[2]);
      arma::vec3 D = get_point(pnum[3]);

      arma::vec3 AB = (B - A);
      AB = AB / arma::norm(AB, 2);
      arma::vec3 AD = D - A;
      AD = AD / arma::norm(AD, 2);
      arma::vec3 normal = arma::cross(AB, AD);
      arma::vec3 PA = A - ref;
      PA = PA / arma::norm(PA, 2);

      if (arma::dot(PA, normal) < 0)
      {
        int aux = pnum[0];
        pnum[0] = pnum[2];
        pnum[2] = aux;
      }

      boundary_elems.push_back(Element(ELEM_QUAD, pnum, marker));
    }
  }

}

inline void Mesh::read_fiber(istream &in, int code,
                             arma::vec3 &f, arma::vec3 &t, arma::vec3 &n)
{
  // Code
  //  0 - orthotropic
  //  1 - transversely isotropic
  //  2 - isotropic
  if (code == 0) {
    in >> f[0] >> f[1] >> f[2]
       >> t[0] >> t[1] >> t[2]
       >> n[0] >> n[1] >> n[2];
  } else {
    t.fill(-1);
    n.fill(-1);
    in >> f[0] >> f[1] >> f[2];
  }
}

arma::vec3 Mesh::get_element_centroid(int elnum) const
{
  Element e = elems[elnum];
  std::vector<int> pnums = e.get_pt_nums();
  std::vector<arma::vec3> p(pnums.size());

  for (uint i = 0; i < pnums.size(); i++)
    p[i] = points[pnums[i]];

  // Triangle only
  arma::vec3 c;
  double x = 0.3 * (points[0](0) + points[1](0) + points[2](0));
  double y = 0.3 * (points[0](1) + points[1](1) + points[2](1));
  c[0] = x;
  c[1] = y;
  return c;
}

Mapping Mesh::get_element_mapping(int elnum) const
{
  int ndim;
  Element e = elems[elnum];
  std::vector<int> pnums = e.get_pt_nums();
  std::vector<arma::vec3> p(pnums.size());

  ElementType et = e.get_type();

  for (uint i = 0; i < pnums.size(); i++) p[i] = points[pnums[i]];

  if (et == ELEM_TETRA || et == ELEM_HEXA)
    ndim = 3;
  else if (et == ELEM_QUAD || et == ELEM_TRIG)
    ndim = 2;
  else
    ndim = 1;

  return Mapping(p, ndim, elems[elnum].get_index());
}

Mapping Mesh::get_element_mapping(int elnum, std::vector<arma::vec3> &x) const
{
  Element e = elems[elnum];
  std::vector<int> pnums = e.get_pt_nums();
  std::vector<arma::vec3> p(pnums.size());

  ElementType et = e.get_type();

  int ndim;

  if (et == ELEM_TETRA || et == ELEM_HEXA)
    ndim = 3;
  else if (et == ELEM_QUAD || et == ELEM_TRIG)
    ndim = 2;
  else
    ndim = 1;

  return Mapping(x, ndim, elems[elnum].get_index());
}

Mapping *Mesh::element_mapping(int elnum, std::vector<arma::vec3> &x) const
{
  Element e = elems[elnum];
  std::vector<int> pnums = e.get_pt_nums();
  std::vector<arma::vec3> p(pnums.size());

  ElementType et = e.get_type();

  int ndim;

  if (et == ELEM_TETRA || et == ELEM_HEXA)
    ndim = 3;
  else if (et == ELEM_QUAD || et == ELEM_TRIG)
    ndim = 2;
  else
    ndim = 1;

  return new Mapping(x, ndim, elems[elnum].get_index());
}

SurfaceMapping Mesh::get_boundary_element_mapping(int elnum) const
{
  int ndim;
  Element e = boundary_elems[elnum];
  std::vector<int> pnums = e.get_pt_nums();
  std::vector<arma::vec3> p(pnums.size());

  for (uint i = 0; i < pnums.size(); i++)
    p[i] = points[pnums[i]];

  if (e.get_type() == ELEM_SEGM)
    ndim = 1;
  else
    ndim = 2;

  return SurfaceMapping(p, ndim, boundary_elems[elnum].get_index());
}

SurfaceMapping Mesh::get_boundary_element_mapping(int elnum,
                                                  std::vector<arma::vec3> &x) const
{
  int ndim;
  Element e = boundary_elems[elnum];
  std::vector<int> pnums = e.get_pt_nums();
  std::vector<arma::vec3> p(pnums.size());

  if (e.get_type() == ELEM_SEGM)
    ndim = 1;
  else
    ndim = 2;

  return SurfaceMapping(x, ndim, boundary_elems[elnum].get_index());
}

void Mesh::reserve_points(const uint num)
{
  try {
    points.reserve(num);
  }
  catch (std::bad_alloc &e) {
    error("std::bad_alloc in Mesh::reserve_points()");
    exit(1);
  }
}

void Mesh::reserve_elements(const uint num)
{
  try {
    elems.reserve(num);
  }
  catch (std::bad_alloc &e) {
    error("std::bad_alloc in Mesh::reserve_elements()");
    exit(1);
  }

}

void Mesh::swap_nodes(int e, int a, int b)
{
  elems[e].swap_nums(a, b);
}

void Mesh::set_tri_positive_volume()
{
  // to implement
}

void Mesh::set_tet_positive_volume()
{
  uint i = 0;
  double hmax = 0.0;
  bool is_tetra = true;
  while (is_tetra && i < elems.size()) {
    if (elems[i].get_type() != ELEM_TETRA)
      is_tetra = false;
    i++;
  }

  if (is_tetra)
    cout << " Is tetrahedral mesh: yes" << endl;
  else
    cout << " Is tetrahedral mesh: no" << endl;
  
  if (is_tetra)
  {
    for (uint i = 0; i < elems.size(); i++)
    {
      std::vector<int> epnums = elems[i].get_pt_nums();
      arma::mat Je(3, 3);
      double px[4] = {0}, py[4] = {0}, pz[4] = {0};
      for (uint j = 0; j < epnums.size(); j++)
      {
        uint k = epnums[j];
	const arma::vec3 p = get_point(k);
        px[j] = p(0);
        py[j] = p(1);
        pz[j] = p(2);
      }
      Je.zeros();

      Je(0,0) = px[1] - px[0];
      Je(0,1) = px[2] - px[0];
      Je(0,2) = px[3] - px[0];

      Je(1,0) = py[1] - py[0];
      Je(1,1) = py[2] - py[0];
      Je(1,2) = py[3] - py[0];

      Je(2,0) = pz[1] - pz[0];
      Je(2,1) = pz[2] - pz[0];
      Je(2,2) = pz[3] - pz[0];

      if (arma::det(Je) < 0){
        swap_nodes(i, 0, 2);
      }
      // find max H in the mesh
      // tet edges
      // 0,1 - 0,2 - 1,2 - 0,3 - 1,3 - 2,3
      arma::vec3 p0 = get_point(epnums[0]);
      arma::vec3 p1 = get_point(epnums[1]);
      arma::vec3 p2 = get_point(epnums[2]);
      arma::vec3 p3 = get_point(epnums[3]);      
      std::vector<double> edge_len;
      edge_len.resize(6);
      edge_len[0] = arma::norm(p0-p1,2);
      edge_len[1] = arma::norm(p0-p2,2);      
      edge_len[2] = arma::norm(p1-p2,2);      
      edge_len[3] = arma::norm(p0-p3,2);      
      edge_len[4] = arma::norm(p1-p3,2);
      edge_len[5] = arma::norm(p2-p3,2);
      double ehval = *std::max_element(edge_len.begin(),edge_len.end());      
      if (ehval > hmax)
        hmax = ehval;      
    }    
  } // if tetra

  cout << " Maximum edge size (h_max): " << hmax << endl;  
}

std::ostream &operator<<(std::ostream &os, const Mesh &mesh)
{
  msg("Mesh information");
  os << " Number of dimensions: " << mesh.get_n_dim();
  os << std::endl;
  os << " Number of nodes: " << mesh.get_n_points();
  os << std::endl;
  os << " Number of elements: " << mesh.get_n_elements();
  os << std::endl;
  os << " Number of boundary elements: " << mesh.get_n_boundary_elements();
  os << std::endl;
  return os;
}
