#include "util.hpp"
#include "pugixml.hpp"

// -----------------------------------------------------------------------------
// Printing utils
// -----------------------------------------------------------------------------

void title(std::string t)
{
  std::cout << std::endl;
  std::cout << "\033[1m" << t << "\033[0m" << std::endl;
  std::cout << std::endl;
}

void msg(std::string m)
{
  std::cout << m << std::endl;
}

std::string line()
{
  std::string l;
  l  = "---------------------------------------";
  l += "---------------------------------------";
  return l;
}

void print_error(std::string file, std::string func, std::string erro)
{
  std::stringstream s;
  s << std::endl << std::endl;
  s << "*** ";
  s << "-------------------------------------";
  s << "-------------------------------------";
  s << std::endl;
  s << "*** E R R O R";
  s << std::endl;
  s << "*** ";
  s << "-------------------------------------";
  s << "-------------------------------------";
  s << std::endl;
  s << "*** " << "File: " << file << std::endl;
  s << "*** " << "Function: " << func << "()" << std::endl;
  s << "*** " << "Error: " << erro << "." << std::endl;
  s << "*** ";
  s << "-------------------------------------";
  s << "-------------------------------------";
  s << std::endl;

  throw std::runtime_error(s.str());
}

// -----------------------------------------------------------------------------
// File utils
// -----------------------------------------------------------------------------

bool file_exists(const char * filename)
{
  ifstream file(filename);
  if(file)
    return true;
  else
    return false;
}

bool file_exists(const std::string & filename)
{
  /**
    The return statement should cast the file object
    to a boolean which is true if the file exists.
    The file is automatically closed at the end of the function scope.
   */
  ifstream file(filename.c_str());
  if(file)
    return true;
  else
    return false;
}

std::string file_extension(const std::string & filename)
{
  std::size_t pos = filename.find_last_of(".");
  std::string extension = "";
  if (pos != std::string::npos)
    extension = filename.substr(pos+1);
  return extension;
}

// -----------------------------------------------------------------------------
// String utils
// -----------------------------------------------------------------------------

void string_split(std::string str, std::string separator, 
		  std::vector<std::string> * results)
{
  std::string::size_type found;
  found = str.find_first_of(separator);

  while(found != std::string::npos)
  {
    if(found > 0)
    {
      results->push_back(str.substr(0,found));
    }
    str = str.substr(found+1);
    found = str.find_first_of(separator);
  }

  if(str.length() > 0)
    results->push_back(str);
}

void trim_white_spaces(std::string & str)
{
  // trim leading spaces
  size_t startpos = str.find_first_not_of(" \t");
  if( string::npos != startpos )
  {
    str = str.substr( startpos );
  }

  // trim trailing spaces
  size_t endpos = str.find_last_not_of(" \t");
  if( string::npos != endpos )
  {
    str = str.substr( 0, endpos+1 );
  }
}

inline void loadbar(uint x, uint n, uint w)
{
    if ( (x != n) && (x % (n/100) != 0) ) return;
 
    float ratio = x/(float)n;
    uint c = ratio * w;
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (uint x=0; x<c; x++) cout << "=";
    for (uint x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

void parse_to_vec3(const std::string & str, arma::vec3 & v)
{
  std::stringstream ss(str);
  ss >> v(0);
  if(ss.peek() == ',')
    ss.ignore();

  ss >> v(1);
  if(ss.peek() == ',')
    ss.ignore();

  ss >> v(2);
  if(ss.peek() == ',')
    ss.ignore();
}

// -----------------------------------------------------------------------------
// XML tools
// -----------------------------------------------------------------------------

void read_xml_section(const std::string & xmlname,
                      const std::string & section,
                      const std::string & subsection,
                      std::map<int,double> & pmap)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(xmlname.c_str());

  if (!result)
  {
    std::cout << "XML parsed with errors, attr value: ["
    << doc.child("mesh").attribute("attr").value() << "]\n";
    std::cout << "Error description: " << result.description() << "\n";
    std::cout << "Error offset: " << result.offset;
    std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
  }

  // element type
  std::string celltype = doc.child("mesh").attribute("celltype").value();

  //
  // reading parameters
  //
  pugi::xml_node nodes = doc.child(section.c_str()).child(subsection.c_str());
  for(pugi::xml_node node = nodes.child("node") ; node ;
      node = node.next_sibling("node"))
  {
    int idx = node.attribute("id").as_int();
    int marker = node.attribute("marker").as_int();
    double val = node.attribute("value").as_double();
    pmap.insert( std::pair<int,double>(marker,val) );
  }

}

// -----------------------------------------------------------------------------
// Arma tools
// -----------------------------------------------------------------------------

void voigtvec(const int nd, const arma::mat & m, arma::mat & v)
{
  //v.zeros();

  if (nd == 2)
  { 
    assert(v.n_rows == 3);
    v(0,0) = m(0,0);
    v(1,0) = m(1,1);
    v(2,0) = m(0,1);
  }
  else if (nd == 3)
  {
    assert(v.n_rows == 6);
    v(0,0) = m(0,0);
    v(1,0) = m(1,1);
    v(2,0) = m(2,2);
    v(3,0) = m(0,1);
    v(4,0) = m(1,2);
    v(5,0) = m(0,2);
  }
}

double contract(const arma::mat & A, const arma::mat & B)
{
  double s = 0.0;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      s = s + A(i,j) * B(i,j);
  return s;
}

double trace(const arma::mat & A)
{
  double tr = 0.0;
  for(int i=0; i<3; i++)
    tr = tr + A(i,i);
  return tr;
}
