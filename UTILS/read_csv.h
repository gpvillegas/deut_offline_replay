#include "parse_utils.h"
#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


vector<double> read_csv(string csv_file="", string col_header=""){

  /*
    brief: reads a comma-separated values (CSV) file and returns a vector
    of requested column by user. 

    assumptions: 
    no white-spaces, 
    comments are denoted by '#', 
    first non-commented line is header (name of columns)
    all rows after 1st row are numeric string values (and will be converted to double)
    
    example of format accepted:
    my_file.csv
    ----------------------------
    # this is my csv file
    #
    # another comment
    header_1, header_2, header_3
    34, 30.1, 3.1e5
    21, 12.2, 1.2
    10, 0.41, 0.005
    ---------------------------

    example of interactive usage:
    -------------------------------
    deuteron@physics:~$ root -l
    root [0] .L read_csv.cpp
    root [1] read_csv("my_file.csv", "header_2")
    (std::vector<double>) { 30.100000, 12.200000, 0.41000000 }
  
    If using from a c++ code, please include this file in header as:
    #include "read_csv.cpp"    
    
    then the fucntion may be called within the code.

  */
  
  
  ifstream myFileStream(csv_file.c_str());

  if(!myFileStream.is_open()){
    cout << Form("File %s failed to open",csv_file.c_str()) << endl;
   
  }

  if( col_header.empty() ) {
    cout << "Empty column header !" << endl;
    cout << csv_file.c_str() << "Please select one of the following column headers --- > " << endl;
  }

  string line;
  vector<string> parsed_header;
  int col_idx;
  string row_str;
  double row_value;
  vector<double> col_vec;
  
  int row_cnt = 0;

  // read the file: line by line
  while(getline(myFileStream, line)) {
    stringstream ss(line);

    // ignore comments
    if(line[0]=='#') continue;

    
    // first non-comment row (headers)
    if (row_cnt==0){

      // parsed vector, whose elemetnts are the column headers
      parsed_header = parse_line(line, ',');

      // loop over each header
      for (unsigned int i=0; i < parsed_header.size(); i++){
	// remove leading/trainling spaces in header names
	parsed_header[i] = trim(parsed_header[i]);

	// user help:  display headers if none detected from user
	if(col_header.empty()){
	  cout << parsed_header[i] << endl;

	  if(i == parsed_header.size() - 1 ) {
	    cout << "reached the end of file " << endl;
	    exit(0);
	  }
	}
	
	if (parsed_header[i] == col_header.c_str() ) {
	  col_idx = i; // store the column index requested by user

	}
      }  
    }


    // loop over each column index
    for (unsigned int i=0; i < parsed_header.size(); i++){

      // get the entire row index
      getline(ss, row_str, ',');

      // require user-requested column index
      if (i==col_idx){

	// convert row string to numeric value for storing in array
	if(row_cnt!=0) {

	  row_value = stod(row_str);

	  // add each numeric row value to the user-requested column vector
	  col_vec.push_back(row_value);

	}

      }
      
    }

    
    // row counter
    row_cnt++;
  }
  
  myFileStream.close();

  // return numerical column data vector 
  return col_vec;
}
