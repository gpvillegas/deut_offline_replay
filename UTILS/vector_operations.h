#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// define vector operations


//______________________________________________
double vsum(std::vector<double> v={}){

  /*brief: returns sum all elements of a vector v */
  
  if (v.empty()) {
    cout << "vector is empty !" << endl;
    exit(0);
  }

  double v_sum;
  v_sum = std::accumulate(v.begin(), v.end(), 0.0);

  return v_sum;
  
}

//_________________________________________________________________
vector<double> vpow(std::vector<double> v={}, signed int k=0){

  // brief: returns vector v raised to power k, where k= . . . -2, -1, 0, 1, 2,... integer
  // k < 0 implies inverse power v^{-k} = 1/v^{k}
  
  if (v.empty()) {
    cout << "vector is empty !" << endl;
    exit(0);
  }
  
  std::vector<double> v_pow;
  
  for(int i=0;i<v.size();++i){

    v_pow.push_back( pow(v[i], k) );
    
  }

  return v_pow;
  
}

//______________________________________________________________________________
vector<double> vmult(std::vector<double> v1={}, std::vector<double> v2={}) {

  // brief: returns vector product of two vectors v1[] * v2[]
  // NOTE: this is standard element by element multiplication (different from "dot" or "cross" product)

  if (v1.empty() || v2.empty() ) {
    cout << "vector v1 or v2 is empty !" << endl;
    exit(0);
  }

  std::vector<double> v_prod;
  
  for(int i=0;i<v1.size();++i){

    v_prod.push_back( v1[i] * v2[i] );
      
  }

  return v_prod;
  
}

//_____________________________________________________________________
vector<double> vscale(std::vector<double> v={}, double k=0) {

  // brief: returns vector v multiplied by scaler k, v[] * k


  if (v.empty() ) {
    cout << "vector v is empty !" << endl;
    exit(0);
  }
  
  std::vector<double> v_scl;

  for(int i=0;i<v.size();++i){
    v_scl.push_back( v[i] * k );
  }

  return v_scl;

}


//___________________________________________________________________________
double vavg(std::vector<double> v){

  // brief: calculate normal average of a vector v,

  if (v.empty()) {
    cout << "vector v is empty !" << endl;
    exit(0);
  }
  
  double v_avg=0;
    
    // calculate standad average
    v_avg = vsum(v) / v.size();

    return v_avg;  
  
}

//______________________________________________________________________________
double vavgw(std::vector<double> v, std::vector<double> v_err, double &err_w){

  // brief: returns the calculate weighted average of a vector v with absolute error v_err, where
  // err_w is the uncertainty of the weighted average passed by reference
  

  std::vector<double> weight;

  double num=0;
  double den=0;
  
  double v_avg=0;
  double v_avg_err=0;


  // define weight as inverse square of error
  weight =  vpow(v_err, -2); // weight = 1/v_err^{2}  
  
  
  // element-wise sum of product:  sum { v_i * weight_i }
  num = vsum( vmult(v, weight) );
  
  // element-wise sum of weights, sum { weight_i }
  den = vsum( weight );
  
  // weighted average: v_weighted = sum { v_i * weight_i } / sum { weight_i }
  v_avg = num / den;
  
  // error of weighted average is: sqrt [ 1. / sum { 1 / v_err^{2} } ] =  sqrt [ 1. / sum { weights } ]
  v_avg_err = sqrt( 1. / vsum( weight) );
  
  err_w =  v_avg_err; // assign the error to the varaible passed by reference by the user

  return v_avg;
  
}
