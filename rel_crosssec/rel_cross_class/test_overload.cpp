#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Gas
{
private:
  vector<double> real;
 
public:
  Gas(vector<double> val) { real  = val; }
  friend Gas operator+(const Gas &c1, const Gas &c2); 
  void print()
  {
    for(int i = 0; i < real.size(); ++i)
      cout<<real.at(i)<<",";
    cout<<endl;
  }
};
 
// note: this function is not a member function nor a friend function!
Gas operator+(const Gas &c1, const Gas &c2)
{
	// use the Cents constructor and operator+(int, int)
        // we don't need direct access to private members here
  vector<double> c3;

  for(int i = 0; i < c1.real.size(); ++i)
    {
      double c3_v = c1.real.at(i) + c2.real.at(i);
      c3.push_back(c3_v);
    }

	return Gas(c3);
}

int main()
{
  vector<double> x1 = {1,2};
  vector<double> x2 = {3,4};
  Gas g1(x1);
  Gas g2(x2);
  g1.print();
  cout<<"Vector 2"<<endl;
  g2.print();
  Gas g3 = g1 + g2;
  cout<<"vector mix"<<endl;
  g3.print();
 
	return 0;
}
