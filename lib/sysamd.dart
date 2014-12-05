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

/// The symamd routine computes an ordering P of a symmetric sparse
/// matrix A such that the Cholesky factorization PAP' = LL' remains
/// sparse.  It is based on a column ordering of a matrix M constructed
/// so that the nonzero pattern of M'M is the same as A.  The matrix A
/// is assumed to be symmetric; only the strictly lower triangular part
/// is accessed.  You must pass your selected memory allocator (usually
/// calloc/free or mxCalloc/mxFree) to symamd, for it to allocate
/// memory for the temporary matrix M.
///
/// @param n number of rows and columns of A. Restriction:  n >= 0.
/// Symamd returns FALSE if n is negative.
/// @param A an integer array of size nnz, where nnz = p [n].
///
/// The row indices of the entries in column c of the matrix are
/// held in A [(p [c]) ... (p [c+1]-1)].  The row indices in a
/// given column c need not be in ascending order, and duplicate
/// row indices may be present.  However, symamd will run faster
/// if the columns are in sorted order with no duplicate entries.
///
/// The matrix is 0-based.  That is, rows are in the range 0 to
/// n-1, and columns are in the range 0 to n-1.  Symamd
/// returns FALSE if any row index is out of range.
///
/// The contents of A are not modified.
/// @param an integer array of size n+1.  On input, it holds the
/// "pointers" for the column form of the matrix A.  Column c of
/// the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
/// entry, p [0], must be zero, and p [c] <= p [c+1] must hold
/// for all c in the range 0 to n-1.  The value p [n] is
/// thus the total number of entries in the pattern of the matrix A.
/// Symamd returns FALSE if these conditions are not met.
///
/// The contents of p are not modified.
/// @param perm On output, if symamd returns TRUE, the array perm holds the
/// permutation P, where perm [0] is the first index in the new
/// ordering, and perm [n-1] is the last.  That is, perm [k] = j
/// means that row and column j of A is the kth column in PAP',
/// where k is in the range 0 to n-1 (perm [0] = j means
/// that row and column j of A are the first row and column in
/// PAP').  The array is used as a workspace during the ordering,
/// which is why it must be of length n+1, not just n.
/// @param knobs parameters (uses defaults if NULL)
/// @param stats Statistics on the ordering, and error status.
/// Symamd returns FALSE if stats is not present.
///
/// stats [0]:  number of dense or empty row and columns ignored
///     (and ordered last in the output permutation
///     perm).  Note that a row/column can become
///     "empty" if it contains only "dense" and/or
///     "empty" columns/rows.
///
/// stats [1]:  (same as stats [0])
///
/// stats [2]:  number of garbage collections performed.
///
/// stats [3]:  status code.  < 0 is an error code.
///       > 1 is a warning or notice.
///
///     0  OK.  Each column of the input matrix contained
///       row indices in increasing order, with no
///       duplicates.
///
///     1  OK, but columns of input matrix were jumbled
///       (unsorted columns or duplicate entries).  Symamd
///       had to do some extra work to sort the matrix
///       first and remove duplicate entries, but it
///       still was able to return a valid permutation
///       (return value of symamd was TRUE).
///
///         stats [4]: highest numbered column that
///           is unsorted or has duplicate
///           entries.
///         stats [5]: last seen duplicate or
///           unsorted row index.
///         stats [6]: number of duplicate or
///           unsorted row indices.
///
///     -1  A is a null pointer
///
///     -2  p is a null pointer
///
///     -3  (unused, see colamd.c)
///
///     -4   n is negative
///
///         stats [4]: n
///
///     -5  number of nonzeros in matrix is negative
///
///         stats [4]: # of nonzeros (p [n]).
///
///     -6  p [0] is nonzero
///
///         stats [4]: p [0]
///
///     -7  (unused)
///
///     -8  a column has a negative number of entries
///
///         stats [4]: column with < 0 entries
///         stats [5]: number of entries in col
///
///     -9  a row index is out of bounds
///
///         stats [4]: column with bad row index
///         stats [5]: bad row index
///         stats [6]: n_row, # of rows of matrx
///
///     -10  out of memory (unable to allocate temporary
///       workspace for M or count arrays using the
///       "allocate" routine passed into symamd).
///
/// Future versions may return more statistics in the stats array.
/// @param allocate pointer to calloc
/// @param release pointer to free
/// @return TRUE if OK, FALSE otherwise
int symamd(int n, List<int> A, List<int> p, List<int> perm, List<double> knobs, List<int> stats) {
  /* === Local variables ============================================== */

  Int32List count; // length of each column of M, and col pointer
  Int32List mark; // mark array for finding duplicate entries
  Int32List M; // row indices of matrix M
  int Mlen; // length of M
  int n_row; // number of rows in M
  int nnz; // number of entries in A
  int i; // row index of A
  int j; // column index of A
  int k; // row index of M
  int mnz; // number of nonzeros in M
  int pp; // index into a column of A
  int last_row; // last row seen in the current column
  int length; // number of nonzeros in a column

  final cknobs = new Float64List(KNOBS); // knobs for colamd
  final default_knobs = new Float64List(KNOBS); // default knobs for colamd

  /* === Check the input arguments ==================================== */

  if (stats == null) {
    _debug0("symamd: stats not present");
    return (_false);
  }
  for (i = 0; i < STATS; i++) {
    stats[i] = 0;
  }
  stats[STATUS] = OK;
  stats[INFO1] = -1;
  stats[INFO2] = -1;

  if (A == null) {
    stats[STATUS] = ERROR_A_NOT_PRESENT;
    _debug0("symamd: A not present");
    return (_false);
  }

  if (p == null) // p is not present
  {
    stats[STATUS] = ERROR_P_NOT_PRESENT;
    _debug0("symamd: p not present");
    return (_false);
  }

  if (n < 0) // n must be >= 0
  {
    stats[STATUS] = ERROR_NCOL_NEGATIVE;
    stats[INFO1] = n;
    _debug0("symamd: n negative $n");
    return (_false);
  }

  nnz = p[n];
  if (nnz < 0) // nnz must be >= 0
  {
    stats[STATUS] = ERROR_NNZ_NEGATIVE;
    stats[INFO1] = nnz;
    _debug0("symamd: number of entries negative $nnz");
    return (_false);
  }

  if (p[0] != 0) {
    stats[STATUS] = ERROR_P0_NONZERO;
    stats[INFO1] = p[0];
    _debug0("symamd: p[0] not zero ${p [0]}");
    return (_false);
  }

  /* === If no knobs, set default knobs =============================== */

  if (knobs == null) {
    setDefaults(default_knobs);
    knobs = default_knobs;
  }

  /* === Allocate count and mark ====================================== */

  try {
    count = new Int32List(n + 1);
  } on OutOfMemoryError catch (e) {
    stats[STATUS] = ERROR_OUT_OF_MEMORY;
    _debug0("symamd: allocate count (size ${n+1}) failed");
    return (_false);
  }

  try {
    mark = new Int32List(n + 1);
  } on OutOfMemoryError catch (e) {
    stats[STATUS] = ERROR_OUT_OF_MEMORY;
    count = null;
    _debug0("symamd: allocate mark (size ${n+1}) failed");
    return (_false);
  }

  /* === Compute column counts of M, check if A is valid ============== */

  stats[INFO3] = 0; // number of duplicate or unsorted row indices

  for (i = 0; i < n; i++) {
    mark[i] = -1;
  }

  for (j = 0; j < n; j++) {
    last_row = -1;

    length = p[j + 1] - p[j];
    if (length < 0) {
      /* column pointers must be non-decreasing */
      stats[STATUS] = ERROR_COL_LENGTH_NEGATIVE;
      stats[INFO1] = j;
      stats[INFO2] = length;
      count = null;
      mark = null;
      _debug0("symamd: col $j negative length $length");
      return (_false);
    }

    for (pp = p[j]; pp < p[j + 1]; pp++) {
      i = A[pp];
      if (i < 0 || i >= n) {
        /* row index i, in column j, is out of bounds */
        stats[STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS;
        stats[INFO1] = j;
        stats[INFO2] = i;
        stats[INFO3] = n;
        count = null;
        mark = null;
        _debug0("symamd: row $i col $j out of bounds");
        return (_false);
      }

      if (i <= last_row || mark[i] == j) {
        /* row index is unsorted or repeated (or both), thus col */
        /* is jumbled.  This is a notice, not an error condition. */
        stats[STATUS] = OK_BUT_JUMBLED;
        stats[INFO1] = j;
        stats[INFO2] = i;
        stats[INFO3]++;
        _debug1("symamd: row $i col $j unsorted/duplicate");
      }

      if (i > j && mark[i] != j) {
        /* row k of M will contain column indices i and j */
        count[i]++;
        count[j]++;
      }

      /* mark the row as having been seen in this column */
      mark[i] = j;

      last_row = i;
    }
  }

  /* === Compute column pointers of M ================================= */

  /* use output permutation, perm, for column pointers of M */
  perm[0] = 0;
  for (j = 1; j <= n; j++) {
    perm[j] = perm[j - 1] + count[j - 1];
  }
  for (j = 0; j < n; j++) {
    count[j] = perm[j];
  }

  /* === Construct M ================================================== */

  mnz = perm[n];
  n_row = mnz ~/ 2;
  Mlen = _recommended(mnz, n_row, n);
  try {
    M = new Int32List(Mlen);
    _debug0("symamd: M is $n_row-by-$n with $mnz entries, Mlen = $Mlen");
  } on OutOfMemoryError catch (e) {
    stats[STATUS] = ERROR_OUT_OF_MEMORY;
    count = null;
    mark = null;
    _debug0("symamd: allocate M (size $Mlen) failed");
    return (_false);
  }

  k = 0;

  if (stats[STATUS] == OK) {
    /* Matrix is OK */
    for (j = 0; j < n; j++) {
      _assert(p[j + 1] - p[j] >= 0);
      for (pp = p[j]; pp < p[j + 1]; pp++) {
        i = A[pp];
        _assert(i >= 0 && i < n);
        if (i > j) {
          /* row k of M contains column indices i and j */
          M[count[i]++] = k;
          M[count[j]++] = k;
          k++;
        }
      }
    }
  } else {
    /* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
    _debug0("symamd: Duplicates in A.");
    for (i = 0; i < n; i++) {
      mark[i] = -1;
    }
    for (j = 0; j < n; j++) {
      _assert(p[j + 1] - p[j] >= 0);
      for (pp = p[j]; pp < p[j + 1]; pp++) {
        i = A[pp];
        _assert(i >= 0 && i < n);
        if (i > j && mark[i] != j) {
          /* row k of M contains column indices i and j */
          M[count[i]++] = k;
          M[count[j]++] = k;
          k++;
          mark[i] = j;
        }
      }
    }
  }

  /* count and mark no longer needed */
  count = null;
  mark = null;
  _assert(k == n_row);

  /* === Adjust the knobs for M ======================================= */

  for (i = 0; i < KNOBS; i++) {
    cknobs[i] = knobs[i];
  }

  /* there are no dense rows in M */
  cknobs[DENSE_ROW] = -1.0;
  cknobs[DENSE_COL] = knobs[DENSE_ROW];

  /* === Order the columns of M ======================================= */

  colamd(n_row, n, Mlen, M, perm, cknobs, stats);

  /* Note that the output permutation is now in perm */

  /* === get the statistics for symamd from colamd ==================== */

  /* a dense column in colamd means a dense row and col in symamd */
  stats[DENSE_ROW] = stats[DENSE_COL];

  /* === Free M ======================================================= */

  M = null;
  _debug0("symamd: done.");
  return (_true);
}
