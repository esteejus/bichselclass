#include <iostream>
#include <functional>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <cmath>

using namespace std;

class gsl_function_pp : public gsl_function
 {
    private:
   std::function<double(double)> _func;
    static double invoke(double x, void *params) {
     return static_cast<gsl_function_pp*>(params)->_func(x);
   }

    public:
 gsl_function_pp(std::function<double(double)> const& func) : _func(func){
      function=&gsl_function_pp::invoke;
      params=this;
    }    
};


int main(){

  //  gsl_function_pp Fp( std::bind(&Class::member_function, &(*this),  std::placeholders::_1) );
  //gsl_function *F = static_cast<gsl_function*>(&Fp);
 
 //Declare them (locally) here
  bool t = true;
   double a1  = 1;
   double a2  = 1.;
 // Declare a lambda function that capture all of them by value or reference
 // no need to write another struct with these 2 parameters + class pointer
   auto ptr2 = [&](double x)->double {
     double f = 0;
     if(t == true)
       f += a1;
     f += a2*x;
     return f;
   };
   auto ptr = [&](double x)->double {return (ptr2(x)*ptr2(x));};
 // Cast to GSL in 3 lines using the wrapper 
 std::function<double(double)> F1(ptr);
 gsl_function_pp F2(F1);
 gsl_function *F = static_cast<gsl_function*>(&F2); 

 double result, error;
 gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
 gsl_integration_qags (F, 0, 2 , 0, 1e-3, 1000,w, &result, &error);    
 gsl_integration_workspace_free (w);

 cout<<"Integral is: "<<result<<endl;

  return 0;
}
