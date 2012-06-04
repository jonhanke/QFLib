

// Some silly vector output operators


// Define the << operator for the vector<T> type
template <class T>
ostream & operator<<(ostream & out, const vector<T> & v) {    
  out << "[ ";
  for(long i=0; i < v.size(); i++)
    out << v[i] << " ";
  out << " ]";
  return out;
}


// Define the << operator for the vector< vector<T> > type
template <class T>
ostream & operator<<(ostream & out, const vector< vector<T> > & vv) {    
  out << "[ " << endl;
  for(long i=0; i < vv.size(); i++)
    out << vv[i] << endl;
  out << " ]" << endl;
  return out;
}





