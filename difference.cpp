
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        std::cerr << "Error! Need at least two files!" << std::endl;
        return 1;
    }
    
    std::cerr << "Opening files: " << argv[0] << " and " << argv[1] << std::endl;
    
    std::ifstream file1(argv[1], std::ifstream::in);
    std::ifstream file2(argv[2], std::ifstream::in);
    
    
    while(file1.good() && file2.good())
    {
        double x1, y1, x2, y2;
        
        file1 >> x1 >> y1;
        file2 >> x2 >> y2;
        
        if(x1 != x2)
        {
            std::cerr << "Error! Files do not match!" << std::endl;
            file1.close();
            file2.close();
            return 1;
        }
        
        std::cout << x1 << " " << (y1 - y2) << std::endl;
    }
    
    file1.close();
    file2.close();
    
    return 0;
}
