part of edu.ufl.cise.colamd;

class Col {

  /// index for A of first row in this column, or DEAD
  /// if column is dead
  int start;

  /// number of rows in this column
  int length;

  int _shared1;
  int _shared2;
  int _shared3;
  int _shared4;

  Col();

  /// number of original columns represented by this
  /// col, if the column is alive
  int get thickness => _shared1;

  void set thickness(int thickness) {
    _shared1 = thickness;
  }

  /// parent in parent tree super-column structure, if
  /// the column is dead
  int get parent => _shared1;

  void set parent(int parent) {
    _shared1 = parent;
  }


  /// the score used to maintain heap, if col is alive
  int get score => _shared2;

  void set score(int score) {
    _shared2 = score;
  }

  /// pivot ordering of this column, if col is dead
  int get order => _shared2;

  void set order(int order) {
    _shared2 = order;
  }


  /// head of a hash bucket, if col is at the head of
  /// a degree list
  int get headhash => _shared3;

  void set headhash(int headhash) {
    _shared3 = headhash;
  }

  /// hash value, if col is not in a degree list
  int get hash => _shared3;

  void set hash(int hash) {
    _shared3 = hash;
  }

  /// previous column in degree list, if col is in a
  /// degree list (but not at the head of a degree list)
  int get prev => _shared3;

  void set prev(int prev) {
    _shared3 = prev;
  }


  /// next column, if col is in a degree list
  int get degree_next => _shared4;

  void set degree_next(int degree_next) {
    _shared4 = degree_next;
  }

  /// next column, if col is in a hash list
  int get hash_next => _shared4;

  void set hash_next(int hash_next) {
    _shared4 = hash_next;
  }

}
