#include <map>
#include <vector>


class PrintTools
{
public:
    PrintTools(){}

    void printArray_int(int* rArray, int rSize, std::string rTitle, int rPID);

//    template<typename T>
    void printVector(std::vector<int> rVector, std::string rTitle, int rPID);

    void printMap_int_int(std::map<int, int> rMap, std::string rTitle, int rPID);


private:

};//class
