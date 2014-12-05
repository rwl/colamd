/// Copyright (c) 1998-2007, Timothy A. Davis, All Rights Reserved.
///
/// This library is free software; you can redistribute it and/or
/// modify it under the terms of the GNU Lesser General Public
/// License as published by the Free Software Foundation; either
/// version 2.1 of the License, or (at your option) any later version.
///
/// This library is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
/// Lesser General Public License for more details.
///
/// You should have received a copy of the GNU Lesser General Public
/// License along with this library; if not, write to the Free Software
/// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
part of edu.ufl.cise.colamd;

/* ========================================================================== */
/* === colamd debugging routines ============================================ */
/* ========================================================================== */

/// At this point, all empty rows and columns are dead.  All live columns
/// are "clean" (containing no dead rows) and simplicial (no supercolumns
/// yet).  Rows may contain dead columns, but all live rows contain at
/// least one live column.
void _debugStructures(int n_row, int n_col, List<Row> rows, List<Col> cols, List<int> A, int n_col2) {
  if (!_ndebug) {
    /* === Local variables ============================================== */

    int i, c, cp, cp_end, len, score, r, rp, rp_end, deg;

    /* === Check A, Row, and Col ======================================== */

    for (c = 0; c < n_col; c++) {
      if (_colIsAlive(cols, c)) {
        len = cols[c].length;
        score = cols[c].score;
        _debug4("initial live col $c $len $score");
        _assert(len > 0);
        _assert(score >= 0);
        _assert(cols[c].thickness == 1);
        cp = cols[c].start;
        cp_end = cp + len;
        while (cp < cp_end) {
          r = A[cp++];
          _assert(_rowIsAlive(rows, r));
        }
      } else {
        i = cols[c].order;
        _assert(i >= n_col2 && i < n_col);
      }
    }

    for (r = 0; r < n_row; r++) {
      if (_rowIsAlive(rows, r)) {
        i = 0;
        len = rows[r].length;
        deg = rows[r].degree;
        _assert(len > 0);
        _assert(deg > 0);
        rp = rows[r].start;
        rp_end = rp + len;
        while (rp < rp_end) {
          c = A[rp++];
          if (_colIsAlive(cols, c)) {
            i++;
          }
        }
        _assert(i > 0);
      }
    }
  }
}


/* ========================================================================== */
/* === debug_deg_lists ====================================================== */
/* ========================================================================== */

/// Prints the contents of the degree lists.  Counts the number of columns
/// in the degree list and compares it to the total it should have.  Also
/// checks the row degrees.
void _debugDegLists(int n_row, int n_col, List<Row> rows, List<Col> cols, List<int> head, int min_score, int should, int max_deg) {
  if (!_ndebug) {
    /* === Local variables ============================================== */

    int deg, col, have, row;

    /* === Check the degree lists ======================================= */

    if (n_col > 10000 && debugLevel <= 0) {
      return;
    }
    have = 0;
    var d = "Degree lists: $min_score\n";
    for (deg = 0; deg <= n_col; deg++) {
      col = head[deg];
      if (col == _empty) {
        continue;
      }
      d += "$deg:";
      while (col != _empty) {
        d += " $col";
        have += cols[col].thickness;
        _assert(_colIsAlive(cols, col));
        col = cols[col].degree_next;
      }
      _debug4(d);
    }
    _debug4("should $should have $have");
    _assert(should == have);

    /* === Check the row degrees ======================================== */

    if (n_row > 10000 && debugLevel <= 0) {
      return;
    }
    for (row = 0; row < n_row; row++) {
      if (_rowIsAlive(rows, row)) {
        _assert(rows[row].degree <= max_deg);
      }
    }
  }
}


/* ========================================================================== */
/* === debug_mark =========================================================== */
/* ========================================================================== */

/// Ensures that the tag_mark is less that the maximum and also ensures that
/// each entry in the mark array is less than the tag mark.
void _debugMark(int n_row, List<Row> rows, int tag_mark, int max_mark) {
  if (!_ndebug) {
    /* === Check the Row marks ========================================== */

    _assert(tag_mark > 0 && tag_mark <= max_mark);
    if (n_row > 10000 && debugLevel <= 0) {
      return;
    }
    for (int r = 0; r < n_row; r++) {
      _assert(rows[r].mark < tag_mark);
    }
  }
}


/* ========================================================================== */
/* === debug_matrix ========================================================= */
/* ========================================================================== */

/// Prints out the contents of the columns and the rows.
void _debugMatrix(int n_row, int n_col, List<Row> rows, List<Col> cols, List<int> A) {
  if (!_ndebug) {
    /* === Local variables ============================================== */

    int r;
    int c;
    int rp;
    int rp_end;
    int cp;
    int cp_end;

    /* === Dump the rows and columns of the matrix ====================== */

    if (debugLevel < 3) {
      return;
    }
    _debug3("DUMP MATRIX:");
    for (r = 0; r < n_row; r++) {
      _debug3("Row $r alive? ${_rowIsAlive (rows, r) ? 1 : 0}");
      if (_rowIsDead(rows, r)) {
        continue;
      }
      _debug3("start ${rows [r].start} length ${rows [r].length} degree ${rows [r].degree}");
      rp = rows[r].start;
      rp_end = rp + rows[r].length;
      while (rp < rp_end) {
        c = A[rp++];
        _debug4("  ${_colIsAlive (cols, c) ? 1 : 0} col $c");
      }
    }

    for (c = 0; c < n_col; c++) {
      _debug3("Col $c alive? ${_colIsAlive (cols, c) ? 1 : 0}");
      if (_colIsDead(cols, c)) {
        continue;
      }
      _debug3("start ${cols [c].start} length ${cols [c].length} shared1 ${cols [c].thickness} shared2 ${cols [c].score}");
      cp = cols[c].start;
      cp_end = cp + cols[c].length;
      while (cp < cp_end) {
        r = A[cp++];
        _debug4("  ${_rowIsAlive (rows, r) ? 1 : 0} row $r");
      }
    }
  }
}
