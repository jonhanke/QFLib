

////////////////////////////
// Define the <= operator //
////////////////////////////
bool boolean_theta::operator<=(const boolean_theta & b) const{

  // Check the precisions and ternary forms are equal
  assert(precision() == b.precision());
  assert(QQ == b.QQ);

  // Check the theta functions are also equal
  for(unsigned long i=0; i<_length(); i++) 
    if (b._theta[i] != (_theta[i] | b._theta[i]))
      return false;

  return true;
}


////////////////////////////
// Define the >= operator //
////////////////////////////
bool boolean_theta::operator>=(const boolean_theta & b) const{

  // Check the precisions and ternary forms are equal
  assert(precision() == b.precision());
  assert(QQ == b.QQ);

  // Check the theta functions are also equal
  for(unsigned long i=0; i<_length(); i++) 
    if (_theta[i] != (_theta[i] | b._theta[i]))
      return false;

  return true;
}



////////////////////////////
// Define the == operator //
////////////////////////////
bool boolean_theta::operator==(const boolean_theta & b) const{

  bool equal_flag = (precision() == b.precision());

  if (equal_flag == false) {
    cout << " The precisions are different... " << endl;
    return equal_flag;  
  }    
  
  
  // Check the ternary forms are also equal
  equal_flag = (equal_flag && (QQ == b.QQ));

  if (equal_flag == false) {
    cout << " The forms are different... " << endl;
    return equal_flag;  
  }    


  // Check the theta functions are also equal
  for(unsigned long i=0; i<_length(); i++) {
    equal_flag = (equal_flag && (_theta[i] == b._theta[i]));
    /*
    if (_theta[i] != b._theta[i]) {
      cout << " They differ at place " << i << endl;
      cout << "     _theta[i] = " << _theta[i] << endl; 
      cout << "   b._theta[i] = " << b._theta[i] << endl; 
      }
    */
  }  


  if (equal_flag == false) {
    cout << " The vectors are different... " << endl;
    return equal_flag;  
  }    
  
  return equal_flag;  
}




/////////////////////////////////////////////////////////////////////
// Copy Constructor -- Needed so that the valarray is resized! =)
////////////////////////////////////////////////////////////////////
void boolean_theta::operator=(const boolean_theta & source) {
   
  // Protect against self-assignment
  if (this != &source) {
    /*
    // Do we really want to copy these??
    THETA_DIR = source.THETA_DIR;
    DIST_DIR = source.DIST_DIR;
    EXEC_DIR = source.EXEC_DIR;
    */

    // Copy the essential attributes
    QQ = source.QQ;
    _theta = source._theta;
    _theta_precision = source._theta_precision;
  }

}
