// set_union example
#include <iostream>     // std::cout
#include <algorithm>    // std::set_union, std::sort
#include <vector>       // std::vector

using namespace std;

int main () {
  vector<double> first = {5,10,15,20,25,50, 55, 30, 35,13.32,13.33};
  vector<double>  second = {50,40,30,20,10,100,5, 13.32};
  std::vector<double> v(first.size() + second.size());                      // 0  0  0  0  0  0  0  0  0  0
  std::vector<double>::iterator it;

  std::sort (first.begin(),first.begin()+first.size());     //  5 10 15 20 25
  std::sort (second.begin(),second.begin()+second.size());   // 10 20 30 40 50

  it=std::set_union (first.begin(), first.begin()+first.size(), second.begin(), second.begin()+second.size(), v.begin());
                                               // 5 10 15 20 25 30 40 50  0  0
  //  v.resize(it-v.begin());                      // 5 10 15 20 25 30 40 50

  std::cout << "The union has " << (v.size()) << " elements:\n";
  for (it=v.begin(); it!=v.end(); ++it)
    std::cout << ' ' << *it;
  std::cout << '\n';

  return 0;
}
