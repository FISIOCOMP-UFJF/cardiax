#include <cstdlib>
#include "../mesh.hpp"

using namespace std;

int main(int argc, char** argv)
{
    cout << "Testing XML reader" << endl;
    
    std::string nome = argv[1];

    Mesh msh;
    msh.read_xml(nome);
    cout << msh << endl;
    cout << "Done" << endl;
    
    return 0;
}

