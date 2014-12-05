part of edu.ufl.cise.colamd;

class Row {

  /// index for A of first col in this row
  int start ;
  /// number of principal columns in this row
  int length ;

  int _shared1 ;
  int _shared2 ;

  Row();

  /// number of principal & non-principal columns in row
  int get degree => _shared1 ;

  void set degree(int degree) {
    _shared1 = degree ;
  }

  /// used as a row pointer in init_rows_cols ()
  int get p => _shared1 ;

  void set p(int p) {
    _shared1 = p ;
  }


  /// for computing set differences and marking dead rows
  int get mark => _shared2 ;

  void set mark(int mark) {
    _shared2 = mark ;
  }

  /// first column in row (used in garbage collection)
  int get first_column => _shared2 ;

  void set first_column(int first_column) {
    _shared2 = first_column ;
  }

}
