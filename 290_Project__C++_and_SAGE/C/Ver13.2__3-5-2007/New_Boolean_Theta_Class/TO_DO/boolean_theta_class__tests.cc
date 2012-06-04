
// Testing Routines
// ----------------

// Make a set of ternary and quaternary forms to use to test things.



// Testing the Read/WriteTernaryThetaBinary Routines
void Test__BooleanTernaryTheta_read_write() {

  cout << endl
       << "==============================================================" 
       << endl << endl;

  cout << "Running Test__BooleanTernaryTheta_read_write()" << endl;
  cout << "----------------------------------------------" << endl;

  // Set the parameters
  long Q[6] = {1, 0, 0, 3, 0, 5};
  unsigned long Precision = 10000;
  char temp_filename[] = "test_1-3-5_10000.bindata";
  
  // Compute and write the first form
  cout << "  Computing the first form." << endl;
  boolean_ternary_theta theta1(Q, Precision);
  theta1.compute();
  cout << "  Finished computing the first form." << endl;
  cout << "  Writing the first form." << endl;  
  theta1.write(temp_filename);
  cout << "  Finished writing the first form." << endl;  


  // Read into the second one
  cout << "  Reading the first form into the second form." << endl;
  boolean_ternary_theta theta2(Q, Precision);
  theta2.read(temp_filename);
  cout << "  Finished reading the first form into the second form." << endl;

  // Remove the testing file
  cout << "    (Removing the temporary file.)" << endl;
  system("rm -f test_1-3-5_10000.bindata");
  cout << "    (Finished removing the temporary file.)" << endl;


  // Then compare them
  if (theta1 == theta2)
    cout << "Passed Test_BooleanTernaryTheta_read_write()." << endl;
  else {
    cout << "Failed Test_BooleanTernaryTheta_read_write()." << endl;
    cout << endl;
    cout << "Diagnostics: " << endl;
    theta1.extended_comparison(theta2);
  }


  cout << endl
       << "==============================================================" 
       << endl << endl;
}


// Test the Magma read against the compute to check the computations are correct
void Test__MagmaRead_vs_compute() {

  cout << endl
       << "==============================================================" 
       << endl << endl;

  cout << "Running Test__MagmaRead_vs_compute()" << endl;
  cout << "------------------------------------" << endl;

  // Set the parameters
  long Q[6] = {1, 0, 0, 3, 0, 5};
  unsigned long Precision = 10000;  // Note : We should also test this for precisions 
                                    // larger and smaller than in the Magma file!


  // Read the Magma theta series file
  cout << "  Reading the MAGMA theta series file." << endl;
  boolean_ternary_theta theta1(Q, Precision);
  theta1.read_magma("../Testing_Files/Magma__1_0_0_3_0_5__10000.txt");  // Relative to THETA_DIR
  cout << "  Finished reading the MAGMA theta series file." << endl;


  // Compute the other theta series directly
  cout << "  Computing the theta series directly." << endl;
  boolean_ternary_theta theta2(Q, Precision);
  theta2.compute();
  cout << "  Finished computing the theta series directly." << endl;


  // Then compare them
  cout << endl;
  if (theta1 == theta2)
    cout << "Passed Test__MagmaRead_vs_compute()." << endl;
  else {
    cout << "Failed Test__MagmaRead_vs_compute()." << endl; 
    cout << endl;

    theta1.extended_comparison(theta2);
  }


  cout << endl
       << "==============================================================" 
       << endl << endl;

}



// Test the compute routine against the compute_distributed routine also
void Test__compute_vs_compute_distributed() {

  cout << endl
       << "==============================================================" 
       << endl << endl;

  cout << "Running Test__compute_vs_compute_distributed()" << endl;
  cout << "----------------------------------------------" << endl;

  // Set the parameters
  long Q[6] = {1, 0, 0, 3, 0, 5};

  /* 
  // This would have 5 slices (-4 --> 0) and no small slices, 
  // but compute_small_slices() always does the first slice!
  unsigned long Precision = 100;
  unsigned long slice_min = 8;
  */

  // /*  // This has 4 slices (-3 --> 0) and one small slice
  unsigned long Precision = 100;
  unsigned long slice_min = 10;
  // */

  /*
  unsigned long Precision = 10000;
  unsigned long slice_min = 100;
  */
  
  /*
  unsigned long Precision = 1000000000;  // 10^9
  unsigned long slice_min = 1000;
  */

  //  cout << "Now1" << endl;

  // Read the Magma theta series file
  boolean_ternary_theta theta1(Q, Precision);
    cout << "Now1.1" << endl;
  theta1.compute_distributed(slice_min);

    cout << "Now2" << endl;

  // Compute the other theta series directly
  boolean_ternary_theta theta2(Q, Precision);
    cout << "Now2.1" << endl;
  theta2.compute();

    cout << "Now3" << endl;

  // Then compare them
  cout << endl;
  if (theta1 == theta2)
    cout << "Passed Test__compute_vs_compute_distributed()." << endl;
  else {
    cout << "Failed Test__compute_vs_compute_distributed()." << endl;
    cout << endl;

    theta1.extended_comparison(theta2);
  }

  cout << endl
       << "==============================================================" 
       << endl << endl;

}




// Test the compute routine against the compute_distributed routine also
void Test__compute_vs_multislice_compute_distributed() {

  cout << endl
       << "==============================================================" 
       << endl << endl;

  cout << "Running Test__compute_vs_multislice_compute_distributed()" << endl;
  cout << "---------------------------------------------------------" << endl;

  // Set the parameters
  long Q[6] = {1, 0, 0, 3, 0, 5};

  /* 
  // This would have 5 slices (-4 --> 0) and no small slices, 
  // but compute_small_slices() always does the first slice!
  unsigned long Precision = 100;
  */

  /*  // This has 4 slices (-3 --> 0) and one small slice
  unsigned long Precision = 100;
  */

  /*
  unsigned long Precision = 10000; // 10^5
  */

  unsigned long Precision = 1000000; // 10^7
  
  /*
  unsigned long Precision = 1000000000;  // 10^9
  */

  //  cout << "Now1" << endl;

  // Read the Magma theta series file
  boolean_ternary_theta theta1(Q, Precision);
    cout << "Now1.1" << endl;
  theta1.compute_distributed_new();

    cout << "Now2" << endl;

  // Compute the other theta series directly
  boolean_ternary_theta theta2(Q, Precision);
    cout << "Now2.1" << endl;
  theta2.compute();

    cout << "Now3" << endl;

  // Then compare them
  cout << endl;
  if (theta1 == theta2)
    cout << "Passed Test__compute_vs_multislice_compute_distributed()." << endl;
  else {
    cout << "Failed Test__compute_vs_multislice_compute_distributed()." << endl;
    cout << endl;

    theta1.extended_comparison(theta2);
  }

  cout << endl
       << "==============================================================" 
       << endl << endl;

}
