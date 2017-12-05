#include <iostream>
#include <vector>

using namespace std;


int main(){

  //  struct combine { vector<vector<double>> data; combine3() : data(3,vector<double>(2)) {}};
  struct combine {
    vector<vector<double>> data;
    vector<vector<double>> data;
    vector<vector<double>> data;
  };
  
  struct combine test;
  vector<vector<double>> entry = { {1,2},{1,3},{1,4},{1,5,6}};

  test.data = entry;
  cout<<test.data.at(3).at(2)<<endl;
  
  return 0;
  
}
