#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

int main()
{

  vector<double> cdf_y = { 0,1e-3,0,0,1e-3, 0, .1, .2,.15, .3, .4, .8, .9 ,.999,.998,.999,1,.98};
  vector<double> cdf_x = { 0,1,2,3,4, 5, 6, 7,8, 9, 10, 11, 12 ,13,14,15,16,17};
  vector<double> cdf_x_mon, cdf_y_mon;//monotomic vectors

      for(int j = 0; j < cdf_y.size(); j++)
	  cout<<cdf_y.at(j)<<"\t";
      cout<<endl;
      for(int j = 0; j < cdf_x.size(); j++)
	  cout<<cdf_x.at(j)<<"\t";


  for( int i = 1; i < cdf_y.size()-2; i++)
    {
      cout<<endl;
      double cur_elm_y = cdf_y.at(i);
      double cur_elm_x = cdf_x.at(i);
      if( cur_elm_y > cdf_y.at(i-1) && cur_elm_y < cdf_y.at(i+1))
	{
	cdf_y_mon.push_back(cur_elm_y);
	cdf_x_mon.push_back(cur_elm_x);
	}
      else
	{
	  if(!(cdf_y.at(i+1) > cur_elm_y && cdf_y.at(i+1) < cdf_y.at(i+2)))
	    {
	      cdf_y_mon.push_back(cur_elm_y);
	      cdf_x_mon.push_back(cur_elm_x);
	      i ++;
	    }

	  //	cdf_y_mon.pop_back();
	  //	cdf_x_mon.pop_back();
	}
      //      prev_elm = cur_elm_y;
      for(int j = 0; j < cdf_y_mon.size(); j++)
	  cout<<cdf_y_mon.at(j)<<" ";
      cout<<endl;
      for(int j = 0; j < cdf_x_mon.size(); j++)
	  cout<<cdf_x_mon.at(j)<<" ";

    }
  /*
  for( int i = 0; i < cdf_y.size(); i++)
    {
      cout<<endl;
      double cur_elm_y = cdf_y.at(i);
      double cur_elm_x = cdf_x.at(i);
      if(cur_elm_y > prev_elm)
	{
	cdf_y_mon.push_back(cur_elm_y);
	cdf_x_mon.push_back(cur_elm_x);
	}
      else
	{
	cdf_y_mon.pop_back();
	cdf_x_mon.pop_back();
	}
      prev_elm = cur_elm_y;
      for(int j = 0; j < cdf_y_mon.size(); j++)
	  cout<<cdf_y_mon.at(j)<<" ";
      cout<<endl;
      for(int j = 0; j < cdf_x_mon.size(); j++)
	  cout<<cdf_x_mon.at(j)<<" ";

    }
*/      
  return 0;
      }
