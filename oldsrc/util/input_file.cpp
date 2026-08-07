#include "input_file.h"

InputFile::InputFile(const char * fname)
{
  std::string temp(fname);
  filename = temp;
  separator= " ";

  is.open(temp.c_str());
  if(!is.is_open())
    std::cerr << "erro " << temp << std::endl;
}

InputFile::InputFile(const std::string & fname)
{
  separator= " ";
  filename = fname;
  is.open(fname.c_str());
}

InputFile::~InputFile()
{
  // do nothing
}

void InputFile::close()
{
  is.close();
}

int InputFile::find_section(const char * name)
{
  std::string str, sname(name);

  is.clear();
  is.seekg(0, std::ios::beg);

  while(is)
  {
    getline(is,str);

    // ignore commentaries
    size_t firstpos = str.find_first_not_of(" ");
    if(firstpos != std::string::npos && str[firstpos] == '#')
    { 
      //cout << "POS = " << firstpos << "\t" 
      //     << "CHAR = " << str[firstpos] << endl;
      //cout<< "COMENTARIO " << str << endl;  
      continue;
    }

    // found
    if (str.find(sname) != std::string::npos)
      return 1; 
  }

  return 0;
}

void InputFile::read_key_value(const std::string & line,
			       std::string & k,
			       std::string & v)
{
  uint found;
  found = line.find_first_of(separator);
  if(found != std::string::npos)
  {
    k = line.substr(0,found);
    v = line.substr(found+1,line.length());
    trim_white_spaces(v);
  }
}

void InputFile::read(const char * key, int & value)
{
  std::string skey(key), k, v;
  std::string str;
  is.clear();
  is.seekg(0, std::ios::beg);
  while(is)
  {
    getline(is,str);

    // ignore commentaries
    size_t firstpos = str.find_first_not_of(" ");
    if(firstpos != std::string::npos && str[firstpos] == '#')
    { 
      //cout << "POS = " << firstpos << "\t" 
      //     << "CHAR = " << str[firstpos] << endl;
      //cout<< "COMENTARIO " << str << endl;  
      continue;
    }

    if (str.find(skey) != std::string::npos)
    {
      read_key_value(str,k,v);
      value = atoi(v.c_str());
    }
  }
}

void InputFile::read(const char * key, double & value)
{
  std::string skey(key), k, v;
  std::string str;

  is.clear();
  is.seekg(0, std::ios::beg);

  while(is)
  {
    getline(is,str);

    // ignore commentaries
    size_t firstpos = str.find_first_not_of(" ");
    if(firstpos != std::string::npos && str[firstpos] == '#')
    { 
      //cout << "POS = " << firstpos << "\t" 
      //     << "CHAR = " << str[firstpos] << endl;
      //cout<< "COMENTARIO " << str << endl;  
      continue;
    }

    if (str.find(skey) != std::string::npos)
    {
      read_key_value(str, k, v);
      value = atof(v.c_str());
    }
  }
}

void InputFile::read(const char * key, std::string & value)
{
  std::string skey(key), k;
  std::string str;

  is.clear();
  is.seekg(0, std::ios::beg);

  while(is)
  {
    getline(is,str);

    // ignore commentaries
    size_t firstpos = str.find_first_not_of(" ");
    if(firstpos != std::string::npos && str[firstpos] == '#')
    {
      //cout << "POS = " << firstpos << "\t" 
      //     << "CHAR = " << str[firstpos] << endl;
      //cout<< "COMENTARIO " << str << endl;      
      continue;
    }

    if (str.find(skey)!=std::string::npos)
      read_key_value(str,k,value);
  }
}

void InputFile::read_array(const char * name, std::vector<double> & array)
{
  uint size;
  std::string str;
  std::vector<std::string> mvec;

  is.clear();
  is.seekg(0, std::ios::beg);

  while(is)
  {
    getline(is,str);

    // ignore commentaries
    size_t firstpos = str.find_first_not_of(" ");
    if(firstpos != std::string::npos && str[firstpos] == '#')
    { 
      //cout << "POS = " << firstpos << "\t" 
      //     << "CHAR = " << str[firstpos] << endl;
      //cout<< "COMENTARIO " << str << endl;  
      continue;
    }

    if (str.find(name)!=std::string::npos)
    {
      string_split(str, " ", &mvec);

      if (mvec.size()==0)
        print_error("input_file.cpp", "read_array", "wrong number of array components");

      size = atoi(mvec[1].c_str());

      for (uint i=2; i<2+size; i++)
      {
        double val = atof(mvec[i].c_str());
        array.push_back(val);
      }
    }
  }
}

void InputFile::read_section(const char * name,
                             std::multimap<int,int> & parmap)
{
  int idx, dir;
  std::string sname(name);
  std::string sname_end = "end_" + sname;
  std::string str, line;
  std::string str_idx, str_dir;

  parmap.clear();

  if(find_section(name))
  {
    do // read until the end of the section
    {
      getline(is,line);

      // break line
      uint found = line.find_first_of(" ");
      if(found != std::string::npos)
      {
        str_idx = line.substr(0,found);
        str_dir = line.substr(found+1,line.length());
      }

      // insert data in map
      idx = atoi(str_idx.c_str());
      dir = atoi(str_dir.c_str());
      parmap.insert( std::pair<int,int>(idx,dir) );

    } while(strcmp(line.c_str(),sname_end.c_str()));
  }

  /*
  is.clear();
  is.seekg(0, std::ios::beg);
  while(is)
  {
    getline(is,str);
    if (str.find(sname)!=std::string::npos)
    {
      // read until the end of the section
      do
      {
        getline(is,line);

        // break line
        uint found = line.find_first_of(" ");
        if(found != std::string::npos)
        {
          str_idx = line.substr(0,found);
          str_dir = line.substr(found+1,line.length());
        }

        // insert data in map
        idx = atoi(str_idx.c_str());
        dir = atoi(str_dir.c_str());
        parmap.insert( std::pair<int,int>(idx,dir) );

      } while(strcmp(line.c_str(),sname_end.c_str()));
    }
  }
  */
}

void InputFile::read_section(const char * name,
                             std::map<int,double> & parmap)
{
  int idx;
  double val;
  std::string sname(name);
  std::string sname_end = "end_" + sname;
  std::string str, line;
  std::string str_idx, str_val;

  parmap.clear();

  if(find_section(name))
  {
    do
    {
      getline(is,line);

      // break line
       uint found = line.find_first_of(" ");
      if(found != std::string::npos)
      {
        str_idx = line.substr(0,found);
        str_val = line.substr(found+1,line.length());       
      }

      // insert data in map
      idx = atoi(str_idx.c_str());
      val = atof(str_val.c_str());
      
      if(strcmp(line.c_str(),sname_end.c_str()) != 0)
        parmap.insert( std::pair<int,double>(idx,val) );

    } while(strcmp(line.c_str(),sname_end.c_str()));

  }
}

void InputFile::read_section(const char * name,
    std::map<int,arma::vec3> & parmap)
{
  int node;
  double vecx, vecy, vecz;
  std::string sname(name);
  std::string sname_end = "end_" + sname;
  std::string str, line, tmp;
  std::vector<std::string> r;

  parmap.clear();

  if(find_section(name))
  {
    do // read until the end of the section
    {
      getline(is,line);
      r.clear();

      string_split(line," ",&r);

      if(r.size() == 4)
      {
        tmp  = r[0];
        node = atoi(tmp.c_str());
        tmp  = r[1];
        vecx = atof(tmp.c_str());
        tmp  = r[2];
        vecy = atof(tmp.c_str());
        tmp  = r[3];
        vecz = atof(tmp.c_str());

        // insert data in map
        arma::vec3 pt;
        pt[0] = vecx; pt[1] = vecy; pt[2] = vecz;
	
        parmap.insert( std::pair<int,arma::vec3>(node,pt) );
      }

    } while(strcmp(line.c_str(),sname_end.c_str()));
  }
}

void InputFile::read_section(const char * name,
                             std::multimap<int,NodalData> & parmap)
{
  int node, dir;
  double val;
  std::string sname(name);
  std::string sname_end = "end_" + sname;
  std::string str, line, tmp;
  std::vector<std::string> r;

  parmap.clear();

  if(find_section(name))
  {
    do // read until the end of the section
    {
      getline(is,line);
      r.clear();
      string_split(line," ",&r);

      if(r.size() == 3)
      {
        tmp  = r[0];
        node = atoi(tmp.c_str());
        tmp  = r[1];
        dir  = atoi(tmp.c_str());
        tmp  = r[2];
        val  = atof(tmp.c_str());

        // insert data in map
        NodalData nd(dir,val);
        parmap.insert( std::pair<int,NodalData>(node,nd) );
      }

    } while(strcmp(line.c_str(),sname_end.c_str()));
  }

}
