/**
 * Copyright (c) 1998-2007, Timothy A. Davis, All Rights Reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 */

library edu.ufl.cise.colamd;

import 'dart:typed_data';
import 'dart:math' as math;

part 'col.dart';
part 'row.dart';

/**
 * COLAMD / SYMAMD
 *
 * colamd:  an approximate minimum degree column ordering algorithm,
 * for LU factorization of symmetric or unsymmetric matrices,
 * QR factorization, least squares, interior point methods for
 * linear programming problems, and other related problems.
 *
 * symamd:  an approximate minimum degree ordering algorithm for Cholesky
 * factorization of symmetric matrices.
 *
 * Purpose:
 *
 * Colamd computes a permutation Q such that the Cholesky factorization of
 * (AQ)'(AQ) has less fill-in and requires fewer floating point operations
 * than A'A.  This also provides a good ordering for sparse partial
 * pivoting methods, P(AQ) = LU, where Q is computed prior to numerical
 * factorization, and P is computed during numerical factorization via
 * conventional partial pivoting with row interchanges.  Colamd is the
 * column ordering method used in SuperLU, part of the ScaLAPACK library.
 * It is also available as built-in function in MATLAB Version 6,
 * available from MathWorks, Inc. (http://www.mathworks.com).  This
 * routine can be used in place of colmmd in MATLAB.
 *
 * Symamd computes a permutation P of a symmetric matrix A such that the
 * Cholesky factorization of PAP' has less fill-in and requires fewer
 * floating point operations than A.  Symamd constructs a matrix M such
 * that M'M has the same nonzero pattern of A, and then orders the columns
 * of M using colmmd.  The column ordering of M is then returned as the
 * row and column ordering P of A.
 *
 * Authors:
 *
 * The authors of the code itself are Stefan I. Larimore and Timothy A.
 * Davis (davis at cise.ufl.edu), University of Florida.  The algorithm was
 * developed in collaboration with John Gilbert, Xerox PARC, and Esmond
 * Ng, Oak Ridge National Laboratory.
 *
 * Acknowledgements:
 *
 * This work was supported by the National Science Foundation, under
 * grants DMS-9504974 and DMS-9803599.
 */

/* ========================================================================== */
/* === COLAMD version ======================================================= */
/* ========================================================================== */

/* COLAMD Version 2.4 and later will include the following definitions.
 * As an example, to test if the version you are using is 2.4 or later:
 *
 *  if (COLAMD_VERSION >= COLAMD_VERSION_CODE (2,4)) ...
 */

const DATE = "Jan 25, 2011";
int versionCode(int main, int sub) => main * 1000 + sub;
const MAIN_VERSION = 2 ;
const SUB_VERSION = 7 ;
const SUBSUB_VERSION = 3 ;
final VERSION = versionCode(MAIN_VERSION, SUB_VERSION) ;

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
const KNOBS = 20 ;

/* number of output statistics.  Only stats [0..6] are currently used. */
const STATS = 20 ;

/* knobs [0] and stats [0]: dense row knob and output statistic. */
const DENSE_ROW = 0 ;

/* knobs [1] and stats [1]: dense column knob and output statistic. */
const DENSE_COL = 1 ;

/* knobs [2]: aggressive absorption */
const AGGRESSIVE = 2 ;

/* stats [2]: memory defragmentation count output statistic */
const DEFRAG_COUNT = 2 ;

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
const STATUS = 3 ;

/* stats [4..6]: error info, or info on jumbled columns */
const int INFO1 = 4 ;
const int INFO2 = 5 ;
const int INFO3 = 6 ;

/* error codes returned in stats [3]: */
const int OK        = (0) ;
const int OK_BUT_JUMBLED      = (1) ;
const int ERROR_A_NOT_PRESENT    = (-1) ;
const int ERROR_P_NOT_PRESENT    = (-2) ;
const int ERROR_NROW_NEGATIVE    = (-3) ;
const int ERROR_NCOL_NEGATIVE    = (-4) ;
const int ERROR_NNZ_NEGATIVE    = (-5) ;
const int ERROR_P0_NONZERO      = (-6) ;
const int ERROR_A_TOO_SMALL    = (-7) ;
const int ERROR_COL_LENGTH_NEGATIVE  = (-8) ;
const int ERROR_ROW_INDEX_OUT_OF_BOUNDS  = (-9) ;
const int ERROR_OUT_OF_MEMORY    = (-10) ;
const int ERROR_INTERNAL    = (-999) ;

/* ========================================================================== */
/* === Scaffolding code definitions  ======================================== */
/* ========================================================================== */

void set debug(bool d) {
  _ndebug = !d;
}
bool _ndebug = true;

void set log(bool l) {
  nprint = !l;
}
bool nprint = true;

/*
Our "scaffolding code" philosophy:  In our opinion, well-written library
code should keep its "debugging" code, and just normally have it turned off
by the compiler so as not to interfere with performance.  This serves
several purposes:

1. assertions act as comments to the reader, telling you what the code
   expects at that point.  All assertions will always be true (unless
   there really is a bug, of course).

2. leaving in the scaffolding code assists anyone who would like to modify
   the code, or understand the algorithm (by reading the debugging output,
   one can get a glimpse into what the code is doing).

3. (gasp!) for actually finding bugs.  This code has been heavily tested
   and "should" be fully functional and bug-free ... but you never know...

     The code will become outrageously slow when debugging is
     enabled.  To control the level of debugging output, set an environment
     variable D to 0 (little), 1 (some), 2, 3, or 4 (lots).  When debugging,
     you should see the following message on the standard output:

       colamd: debug version, D = 1 (THIS WILL BE SLOW!)

     or a similar message for symamd.  If you don't, then debugging has not
     been enabled.
*/

/* ========================================================================== */
/* === Definitions ========================================================== */
/* ========================================================================== */

//double sqrt (double a)
//{
//  return math.sqrt (a) ;
//}

num _max(num a, num b) => a > b ? a : b;

num _min(num a, num b) => a < b ? a : b;

int _denseDegree (double alpha, int n) {
  return _max(16.0, (alpha) * math.sqrt(n.toDouble())).toInt();
}

int _onesComplement(int r) => -(r)-1;

const int _true = 1;
const int _false = 0;

const int _empty = -1;

/* Row and column status */
const int _alive = 0;
const int _dead = -1;

/* Column status */
const int deadPrincipal = -1;
const int deadNonPrincipal = -2;

/* Row and column status update and checking. */
bool _rowIsDead(List<Row> Row, int r)
{
  return _rowIsMarkedDead (Row [r].mark) ;
}
bool _rowIsMarkedDead(int row_mark)
{
  return (row_mark < _alive) ;
}
bool _rowIsAlive(List<Row> Row, int r)
{
  return (Row [r].mark >= _alive) ;
}
bool _colIsDead(List<Col> Col, int c)
{
  return (Col [c].start < _alive) ;
}
bool _colIsAlive(List<Col> Col, int c)
{
  return (Col [c].start >= _alive) ;
}
bool _colIsDeadPrincipal(List<Col> Col, int c)
{
  return (Col [c].start == deadPrincipal) ;
}
void _killRow(List<Row> Row, int r)
{
  Row [r].mark = _dead ;
}
void _killPrincipalCol(List<Col> Col, int c)
{
  Col [c].start = deadPrincipal ;
}
void _killNonPrincipalCol(List<Col> Col, int c)
{
  Col [c].start = deadNonPrincipal ;
}

/* ========================================================================== */
/* === Colamd reporting mechanism =========================================== */
/* ========================================================================== */

/// In Dart, matrices are 0-based and indices are reported as such
/// in *_report.
int _index(int i) => i;

/// All output goes through PRINTF.
void _print(String format) {
  if (!nprint) {
    print(format);
  }
}

/** debug print level */
int debugLevel = 0 ;

void _debug0(String format)
{
  if (!_ndebug)
  {
    _print (format);
  }
}

void _debug1(String format) {
  if (!_ndebug) {
    if (debugLevel >= 1) _print (format);
  }
}

void _debug2(String format)
{
  if (!_ndebug)
  {
    if (debugLevel >= 2) _print (format);
  }
}

void _debug3(String format)
{
  if (!_ndebug)
  {
    if (debugLevel >= 3) _print (format);
  }
}

void _debug4(String format) {
  if (!_ndebug)
  {
    if (debugLevel >= 4) _print (format);
  }
}

void _assert(bool a)
{
  if (!_ndebug)
  {
    assert(a) ;
  }
}

void _assertInt(int a)
{
  _assert (a != 0) ;
}

/* ========================================================================== */
/* === USER-CALLABLE ROUTINES: ============================================== */
/* ========================================================================== */

/* ========================================================================== */
/* === colamd_recommended =================================================== */
/* ========================================================================== */

/*
    The colamd_recommended routine returns the suggested size for Alen.  This
    value has been determined to provide good balance between the number of
    garbage collections and the memory requirements for colamd.  If any
    argument is negative, or if integer overflow occurs, a 0 is returned as an
    error condition.  2*nnz space is required for the row and column
    indices of the matrix. COLAMD_C (n_col) + COLAMD_R (n_row) space is
    required for the Col and Row arrays, respectively, which are internal to
    colamd (roughly 6*n_col + 4*n_row).  An additional n_col space is the
    minimal amount of "elbow room", and nnz/5 more space is recommended for
    run time efficiency.

    Alen is approximately 2.2*nnz + 7*n_col + 4*n_row + 10.

    This function is not needed when using symamd.
*/

/// add two values of type int, and check for integer overflow
int _add(int a, int b, Int32List ok)
{
    ok[0] = (ok[0] != 0) && ((a + b) >= _max(a,b)) ? 1 : 0;
    return ok[0] != 0 ? a + b : 0;
}

/// compute a*k where k is a small integer, and check for integer overflow
int _mult (int a, int k, Int32List ok)
{
    int i, s = 0 ;
    for (i = 0 ; i < k ; i++)
    {
  s = _add (s, a, ok) ;
    }
    return (s) ;
}

/// size of the Col and Row structures
int _c(int n_col, Int32List ok)
{
//    return ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (int))) ;
  return _add (n_col, 1, ok) ;
}

int _r(int n_row, Int32List ok)
{
//      return ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (int))) ;
    return _add (n_row, 1, ok) ;
}

/// @param nnz number of nonzeros in A. This must
/// be the same value as p [n_col] in the call to
/// colamd - otherwise you will get a wrong value
/// of the recommended memory to use.
/// @param n_row number of rows in A
/// @param n_col number of columns in A
/// @return recommended value of Alen. Returns 0
/// if any input argument is negative.  The use of this routine
/// is optional.  Not needed for symamd, which dynamically allocates
/// its own memory.
int _recommended (int nnz, int n_row, int n_col)
{
  int s ;
  Int32List ok = new Int32List.fromList([ _true ]) ;
  if (nnz < 0 || n_row < 0 || n_col < 0)
  {
    return (0) ;
  }
  s = _mult (nnz, 2, ok) ;      /* 2*nnz */
//    c = COLAMD_C (n_col, ok) ;      /* size of column structures */
//    r = COLAMD_R (n_row, ok) ;      /* size of row structures */
//    s = t_add (s, c, ok) ;
//    s = t_add (s, r, ok) ;
  s = _add (s, n_col, ok) ;      /* elbow room */
  s = _add (s, nnz~/5, ok) ;      /* elbow room */
  ok[0] = 1;//(s < Int_MAX) ? 1 : 0;
  return (ok[0] != 0 ? s : 0) ;
}


/* ========================================================================== */
/* === colamd_set_defaults ================================================== */
/* ========================================================================== */

/*
    The colamd_set_defaults routine sets the default values of the user-
    controllable parameters for colamd and symamd:

  Colamd: rows with more than max (16, knobs [0] * sqrt (n_col))
  entries are removed prior to ordering.  Columns with more than
  max (16, knobs [1] * sqrt (MIN (n_row,n_col))) entries are removed
  prior to ordering, and placed last in the output column ordering.

  Symamd: Rows and columns with more than max (16, knobs [0] * sqrt (n))
  entries are removed prior to ordering, and placed last in the
  output ordering.

  knobs [0]  dense row control

  knobs [1]  dense column control

  knobs [2]  if nonzero, do aggresive absorption

  knobs [3..19]  unused, but future versions might use this

*/

/**
 * knobs [0] and knobs [1] control dense row and col detection:
 *
 * Colamd: rows with more than
 * max (16, knobs [DENSE_ROW] * sqrt (n_col))
 * entries are removed prior to ordering.  Columns with more than
 * max (16, knobs [DENSE_COL] * sqrt (MIN (n_row,n_col)))
 * entries are removed prior to
 * ordering, and placed last in the output column ordering.
 *
 * Symamd: uses only knobs [DENSE_ROW], which is knobs [0].
 * Rows and columns with more than
 * max (16, knobs [DENSE_ROW] * sqrt (n))
 * entries are removed prior to ordering, and placed last in the
 * output ordering.
 *
 * COLAMD_DENSE_ROW and COLAMD_DENSE_COL are defined as 0 and 1,
 * respectively, in colamd.h.  Default values of these two knobs
 * are both 10.  Currently, only knobs [0] and knobs [1] are
 * used, but future versions may use more knobs.  If so, they will
 * be properly set to their defaults by the future version of
 * colamd_set_defaults, so that the code that calls colamd will
 * not need to change, assuming that you either use
 * colamd_set_defaults, or pass a (double *) NULL pointer as the
 * knobs array to colamd or symamd.
 *
 * knobs [2]: aggressive absorption
 *
 * knobs [AGGRESSIVE] controls whether or not to do
 * aggressive absorption during the ordering.  Default is TRUE.
 *
 * @param knobs knob array
 */
void setDefaults (Float64List knobs)
{
  /* === Local variables ============================================== */

  int i ;

  if (knobs == null || knobs.length == 0)
  {
    return ;      // no knobs to initialize
  }
  for (i = 0 ; i < KNOBS ; i++)
  {
    knobs [i] = 0.0 ;
  }
  knobs [DENSE_ROW] = 10.0 ;
  knobs [DENSE_COL] = 10.0 ;
  knobs [AGGRESSIVE] = _true.toDouble() ;  // default: do aggressive absorption
}

/* ========================================================================== */
/* === symamd =============================================================== */
/* ========================================================================== */

/**
 * The symamd routine computes an ordering P of a symmetric sparse
 * matrix A such that the Cholesky factorization PAP' = LL' remains
 * sparse.  It is based on a column ordering of a matrix M constructed
 * so that the nonzero pattern of M'M is the same as A.  The matrix A
 * is assumed to be symmetric; only the strictly lower triangular part
 * is accessed.  You must pass your selected memory allocator (usually
 * calloc/free or mxCalloc/mxFree) to symamd, for it to allocate
 * memory for the temporary matrix M.
 *
 * @param n number of rows and columns of A. Restriction:  n >= 0.
 * Symamd returns FALSE if n is negative.
 * @param A an integer array of size nnz, where nnz = p [n].
 *
 * The row indices of the entries in column c of the matrix are
 * held in A [(p [c]) ... (p [c+1]-1)].  The row indices in a
 * given column c need not be in ascending order, and duplicate
 * row indices may be present.  However, symamd will run faster
 * if the columns are in sorted order with no duplicate entries.
 *
 * The matrix is 0-based.  That is, rows are in the range 0 to
 * n-1, and columns are in the range 0 to n-1.  Symamd
 * returns FALSE if any row index is out of range.
 *
 * The contents of A are not modified.
 * @param an integer array of size n+1.  On input, it holds the
 * "pointers" for the column form of the matrix A.  Column c of
 * the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
 * entry, p [0], must be zero, and p [c] <= p [c+1] must hold
 * for all c in the range 0 to n-1.  The value p [n] is
 * thus the total number of entries in the pattern of the matrix A.
 * Symamd returns FALSE if these conditions are not met.
 *
 * The contents of p are not modified.
 * @param perm On output, if symamd returns TRUE, the array perm holds the
 * permutation P, where perm [0] is the first index in the new
 * ordering, and perm [n-1] is the last.  That is, perm [k] = j
 * means that row and column j of A is the kth column in PAP',
 * where k is in the range 0 to n-1 (perm [0] = j means
 * that row and column j of A are the first row and column in
 * PAP').  The array is used as a workspace during the ordering,
 * which is why it must be of length n+1, not just n.
 * @param knobs parameters (uses defaults if NULL)
 * @param stats Statistics on the ordering, and error status.
 * Symamd returns FALSE if stats is not present.
 *
 * stats [0]:  number of dense or empty row and columns ignored
 *     (and ordered last in the output permutation
 *     perm).  Note that a row/column can become
 *     "empty" if it contains only "dense" and/or
 *     "empty" columns/rows.
 *
 * stats [1]:  (same as stats [0])
 *
 * stats [2]:  number of garbage collections performed.
 *
 * stats [3]:  status code.  < 0 is an error code.
 *       > 1 is a warning or notice.
 *
 *     0  OK.  Each column of the input matrix contained
 *       row indices in increasing order, with no
 *       duplicates.
 *
 *     1  OK, but columns of input matrix were jumbled
 *       (unsorted columns or duplicate entries).  Symamd
 *       had to do some extra work to sort the matrix
 *       first and remove duplicate entries, but it
 *       still was able to return a valid permutation
 *       (return value of symamd was TRUE).
 *
 *         stats [4]: highest numbered column that
 *           is unsorted or has duplicate
 *           entries.
 *         stats [5]: last seen duplicate or
 *           unsorted row index.
 *         stats [6]: number of duplicate or
 *           unsorted row indices.
 *
 *     -1  A is a null pointer
 *
 *     -2  p is a null pointer
 *
 *     -3  (unused, see colamd.c)
 *
 *     -4   n is negative
 *
 *         stats [4]: n
 *
 *     -5  number of nonzeros in matrix is negative
 *
 *         stats [4]: # of nonzeros (p [n]).
 *
 *     -6  p [0] is nonzero
 *
 *         stats [4]: p [0]
 *
 *     -7  (unused)
 *
 *     -8  a column has a negative number of entries
 *
 *         stats [4]: column with < 0 entries
 *         stats [5]: number of entries in col
 *
 *     -9  a row index is out of bounds
 *
 *         stats [4]: column with bad row index
 *         stats [5]: bad row index
 *         stats [6]: n_row, # of rows of matrx
 *
 *     -10  out of memory (unable to allocate temporary
 *       workspace for M or count arrays using the
 *       "allocate" routine passed into symamd).
 *
 * Future versions may return more statistics in the stats array.
 * @param allocate pointer to calloc
 * @param release pointer to free
 * @return TRUE if OK, FALSE otherwise
 */
int symamd (int n, List<int> A, List<int> p, List<int> perm,
    List<double> knobs, List<int> stats)
{
  /* === Local variables ============================================== */

  Int32List count ;   // length of each column of M, and col pointer
  Int32List mark ;    // mark array for finding duplicate entries
  Int32List M ;       // row indices of matrix M
  int Mlen ;      // length of M
  int n_row ;     // number of rows in M
  int nnz ;       // number of entries in A
  int i ;         // row index of A
  int j ;         // column index of A
  int k ;         // row index of M
  int mnz ;       // number of nonzeros in M
  int pp ;        // index into a column of A
  int last_row ;  // last row seen in the current column
  int length ;    // number of nonzeros in a column

  final cknobs = new Float64List(KNOBS) ;         // knobs for colamd
  final default_knobs = new Float64List(KNOBS) ;  // default knobs for colamd

  if (!_ndebug)
  {
//      colamd_get_debug ("symamd") ;
  }

  /* === Check the input arguments ==================================== */

  if (stats == null)
  {
    _debug0 ("symamd: stats not present\n") ;
    return (_false) ;
  }
  for (i = 0 ; i < STATS ; i++)
  {
    stats [i] = 0 ;
  }
  stats [STATUS] = OK ;
  stats [INFO1] = -1 ;
  stats [INFO2] = -1 ;

  if (A == null)
  {
    stats [STATUS] = ERROR_A_NOT_PRESENT ;
    _debug0 ("symamd: A not present\n") ;
    return (_false) ;
  }

  if (p == null)    /* p is not present */
  {
    stats [STATUS] = ERROR_P_NOT_PRESENT ;
    _debug0 ("symamd: p not present\n") ;
    return (_false) ;
  }

  if (n < 0)    /* n must be >= 0 */
  {
    stats [STATUS] = ERROR_NCOL_NEGATIVE ;
    stats [INFO1] = n ;
    _debug0 ("symamd: n negative $n\n") ;
    return (_false) ;
  }

  nnz = p [n] ;
  if (nnz < 0)  /* nnz must be >= 0 */
  {
    stats [STATUS] = ERROR_NNZ_NEGATIVE ;
    stats [INFO1] = nnz ;
    _debug0 ("symamd: number of entries negative $nnz\n") ;
    return (_false) ;
  }

  if (p [0] != 0)
  {
    stats [STATUS] = ERROR_P0_NONZERO ;
    stats [INFO1] = p [0] ;
    _debug0 ("symamd: p[0] not zero ${p [0]}\n") ;
    return (_false) ;
  }

  /* === If no knobs, set default knobs =============================== */

  if (knobs == null)
  {
    setDefaults (default_knobs) ;
    knobs = default_knobs ;
  }

  /* === Allocate count and mark ====================================== */

  try
  {
    count = new Int32List(n+1) ;
  } on OutOfMemoryError catch (e) {
    stats [STATUS] = ERROR_OUT_OF_MEMORY ;
    _debug0 ("symamd: allocate count (size ${n+1}) failed\n") ;
    return (_false) ;
  }

  try
  {
    mark = new Int32List(n+1) ;
  } on OutOfMemoryError catch (e) {
    stats [STATUS] = ERROR_OUT_OF_MEMORY ;
    count = null ;
    _debug0 ("symamd: allocate mark (size ${n+1}) failed\n") ;
    return (_false) ;
  }

  /* === Compute column counts of M, check if A is valid ============== */

  stats [INFO3] = 0 ;  // number of duplicate or unsorted row indices

  for (i = 0 ; i < n ; i++)
  {
    mark [i] = -1 ;
  }

  for (j = 0 ; j < n ; j++)
  {
    last_row = -1 ;

    length = p [j+1] - p [j] ;
    if (length < 0)
    {
      /* column pointers must be non-decreasing */
      stats [STATUS] = ERROR_COL_LENGTH_NEGATIVE ;
      stats [INFO1] = j ;
      stats [INFO2] = length ;
      count = null ;
      mark = null ;
      _debug0 ("symamd: col $j negative length $length\n") ;
      return (_false) ;
    }

    for (pp = p [j] ; pp < p [j+1] ; pp++)
    {
      i = A [pp] ;
      if (i < 0 || i >= n)
      {
        /* row index i, in column j, is out of bounds */
        stats [STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS ;
        stats [INFO1] = j ;
        stats [INFO2] = i ;
        stats [INFO3] = n ;
        count = null ;
        mark = null ;
        _debug0 ("symamd: row $i col $j out of bounds\n") ;
        return (_false) ;
      }

      if (i <= last_row || mark [i] == j)
      {
        /* row index is unsorted or repeated (or both), thus col */
        /* is jumbled.  This is a notice, not an error condition. */
        stats [STATUS] = OK_BUT_JUMBLED ;
        stats [INFO1] = j ;
        stats [INFO2] = i ;
        stats [INFO3] ++ ;
        _debug1 ("symamd: row $i col $j unsorted/duplicate\n") ;
      }

      if (i > j && mark [i] != j)
      {
        /* row k of M will contain column indices i and j */
        count [i]++ ;
        count [j]++ ;
      }

      /* mark the row as having been seen in this column */
      mark [i] = j ;

      last_row = i ;
    }
  }

  /* === Compute column pointers of M ================================= */

  /* use output permutation, perm, for column pointers of M */
  perm [0] = 0 ;
  for (j = 1 ; j <= n ; j++)
  {
    perm [j] = perm [j-1] + count [j-1] ;
  }
  for (j = 0 ; j < n ; j++)
  {
    count [j] = perm [j] ;
  }

  /* === Construct M ================================================== */

  mnz = perm [n] ;
  n_row = mnz ~/ 2 ;
  Mlen = _recommended (mnz, n_row, n) ;
  try
  {
    M = new Int32List(Mlen) ;
    _debug0 ("symamd: M is $n_row-by-$n with $mnz entries, Mlen = $Mlen\n") ;
  } on OutOfMemoryError catch (e) {
    stats [STATUS] = ERROR_OUT_OF_MEMORY ;
    count = null ;
    mark = null;
    _debug0 ("symamd: allocate M (size $Mlen) failed\n") ;
    return (_false) ;
  }

  k = 0 ;

  if (stats [STATUS] == OK)
  {
    /* Matrix is OK */
    for (j = 0 ; j < n ; j++)
    {
      _assert (p [j+1] - p [j] >= 0) ;
      for (pp = p [j] ; pp < p [j+1] ; pp++)
      {
        i = A [pp] ;
        _assert (i >= 0 && i < n) ;
        if (i > j)
        {
          /* row k of M contains column indices i and j */
          M [count [i]++] = k ;
          M [count [j]++] = k ;
          k++ ;
        }
      }
    }
  }
  else
  {
    /* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
    _debug0 ("symamd: Duplicates in A.\n") ;
    for (i = 0 ; i < n ; i++)
    {
      mark [i] = -1 ;
    }
    for (j = 0 ; j < n ; j++)
    {
      _assert (p [j+1] - p [j] >= 0) ;
      for (pp = p [j] ; pp < p [j+1] ; pp++)
      {
        i = A [pp] ;
        _assert (i >= 0 && i < n) ;
        if (i > j && mark [i] != j)
        {
          /* row k of M contains column indices i and j */
          M [count [i]++] = k ;
          M [count [j]++] = k ;
          k++ ;
          mark [i] = j ;
        }
      }
    }
  }

  /* count and mark no longer needed */
  count = null ;
  mark = null ;
  _assert (k == n_row) ;

  /* === Adjust the knobs for M ======================================= */

  for (i = 0 ; i < KNOBS ; i++)
  {
    cknobs [i] = knobs [i] ;
  }

  /* there are no dense rows in M */
  cknobs [DENSE_ROW] = -1.0 ;
  cknobs [DENSE_COL] = knobs [DENSE_ROW] ;

  /* === Order the columns of M ======================================= */

  colamd (n_row, n, Mlen, M, perm, cknobs, stats) ;

  /* Note that the output permutation is now in perm */

  /* === get the statistics for symamd from colamd ==================== */

  /* a dense column in colamd means a dense row and col in symamd */
  stats [DENSE_ROW] = stats [DENSE_COL] ;

  /* === Free M ======================================================= */

  M = null ;
  _debug0 ("symamd: done.\n") ;
  return (_true) ;
}

/* ========================================================================== */
/* === colamd =============================================================== */
/* ========================================================================== */

/*
The colamd routine computes a column ordering Q of a sparse matrix
A such that the LU factorization P(AQ) = LU remains sparse, where P is
selected via partial pivoting.   The routine can also be viewed as
providing a permutation Q such that the Cholesky factorization
(AQ)'(AQ) = LL' remains sparse.
*/

/// Computes a column ordering (Q) of A such that P(AQ)=LU or
/// (AQ)'AQ=LL' have less fill-in and require fewer floating point
/// operations than factorizing the unpermuted matrix A or A'A,
/// respectively.
///
/// @param n_row number of rows in A. Restriction:  n_row >= 0.
/// Colamd returns FALSE if n_row is negative.
/// @param n_col number of columns in A. Restriction:  n_col >= 0.
/// Colamd returns FALSE if n_col is negative.
/// @param Alen length of A. Restriction (see note):
/// Alen >= 2*nnz + 6*(n_col+1) + 4*(n_row+1) + n_col
/// Colamd returns FALSE if these conditions are not met.
///
/// Note:  this restriction makes an modest assumption regarding
/// the size of the two typedef's structures in colamd.h.
/// We do, however, guarantee that
///
///     Alen >= colamd_recommended (nnz, n_row, n_col)
///
/// will be sufficient.  Note: the macro version does not check
/// for integer overflow, and thus is not recommended.  Use
/// the colamd_recommended routine instead.
/// @param A row indices of A.
///
/// A is an integer array of size Alen.  Alen must be at least as
/// large as the bare minimum value given above, but this is very
/// low, and can result in excessive run time.  For best
/// performance, we recommend that Alen be greater than or equal to
/// colamd_recommended (nnz, n_row, n_col), which adds
/// nnz/5 to the bare minimum value given above.
///
/// On input, the row indices of the entries in column c of the
/// matrix are held in A [(p [c]) ... (p [c+1]-1)].  The row indices
/// in a given column c need not be in ascending order, and
/// duplicate row indices may be be present.  However, colamd will
/// work a little faster if both of these conditions are met
/// (Colamd puts the matrix into this format, if it finds that the
/// the conditions are not met).
///
/// The matrix is 0-based.  That is, rows are in the range 0 to
/// n_row-1, and columns are in the range 0 to n_col-1.  Colamd
/// returns FALSE if any row index is out of range.
///
/// The contents of A are modified during ordering, and are
/// undefined on output.
/// @param p pointers to columns in A.
///
/// p is an integer array of size n_col+1.  On input, it holds the
/// "pointers" for the column form of the matrix A.  Column c of
/// the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
/// entry, p [0], must be zero, and p [c] <= p [c+1] must hold
/// for all c in the range 0 to n_col-1.  The value p [n_col] is
/// thus the total number of entries in the pattern of the matrix A.
/// Colamd returns FALSE if these conditions are not met.
///
/// On output, if colamd returns TRUE, the array p holds the column
/// permutation (Q, for P(AQ)=LU or (AQ)'(AQ)=LL'), where p [0] is
/// the first column index in the new ordering, and p [n_col-1] is
/// the last.  That is, p [k] = j means that column j of A is the
/// kth pivot column, in AQ, where k is in the range 0 to n_col-1
/// (p [0] = j means that column j of A is the first column in AQ).
///
/// If colamd returns FALSE, then no permutation is returned, and
/// p is undefined on output.
/// @param knobs parameters (uses defaults if NULL)
/// @param stats output statistics and error codes.
///
/// Statistics on the ordering, and error status.
/// Colamd returns FALSE if stats is not present.
///
/// stats [0]:  number of dense or empty rows ignored.
///
/// stats [1]:  number of dense or empty columns ignored (and
///     ordered last in the output permutation p)
///     Note that a row can become "empty" if it
///     contains only "dense" and/or "empty" columns,
///     and similarly a column can become "empty" if it
///     only contains "dense" and/or "empty" rows.
///
/// stats [2]:  number of garbage collections performed.
///     This can be excessively high if Alen is close
///     to the minimum required value.
///
/// stats [3]:  status code.  < 0 is an error code.
///       > 1 is a warning or notice.
///
///     0  OK.  Each column of the input matrix contained
///       row indices in increasing order, with no
///       duplicates.
///
///     1  OK, but columns of input matrix were jumbled
///       (unsorted columns or duplicate entries).  Colamd
///       had to do some extra work to sort the matrix
///       first and remove duplicate entries, but it
///       still was able to return a valid permutation
///       (return value of colamd was TRUE).
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
///     -3   n_row is negative
///
///         stats [4]: n_row
///
///     -4  n_col is negative
///
///         stats [4]: n_col
///
///     -5  number of nonzeros in matrix is negative
///
///         stats [4]: number of nonzeros, p [n_col]
///
///     -6  p [0] is nonzero
///
///         stats [4]: p [0]
///
///     -7  A is too small
///
///         stats [4]: required size
///         stats [5]: actual size (Alen)
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
///     -10  (unused; see symamd.c)
///
///     -999  (unused; see symamd.c)
///
/// Future versions may return more statistics in the stats array.
/// @return TRUE if successful, FALSE otherwise
int colamd (int n_row, int n_col, int Alen, List<int> A,
            List<int> p, List<double> knobs, List<int> stats)
{
  /* === Local variables ============================================== */

  //int i ;      /* loop index */
  int nnz ;      /* nonzeros in A */
  int Row_size ;    /* size of Row [], in integers */
  int Col_size ;    /* size of Col [], in integers */
  int need ;    /* minimum required length of A */
  List<Row> rows ;    /* pointer into A of Row [0..n_row] array */
  List<Col> cols ;    /* pointer into A of Col [0..n_col] array */
  Int32List n_col2 = new Int32List(1) ;    /* number of non-dense, non-empty columns */
  Int32List n_row2 = new Int32List(1) ;    /* number of non-dense, non-empty rows */
  int ngarbage ;    /* number of garbage collections performed */
  Int32List max_deg = new Int32List(1) ;    /* maximum row degree */
  Float64List default_knobs = new Float64List(KNOBS) ;  /* default knobs array */
  int aggressive ;    /* do aggressive absorption */
  Int32List ok ;

  if (!_ndebug)
  {
//      colamd_get_debug ("colamd") ;
  }

  /* === Check the input arguments ==================================== */

  if (stats == null)
  {
    _debug0 ("colamd: stats not present\n") ;
    return (_false) ;
  }
  for (int i = 0 ; i < STATS ; i++)
  {
    stats [i] = 0 ;
  }
  stats [STATUS] = OK ;
  stats [INFO1] = -1 ;
  stats [INFO2] = -1 ;

  if (A == null)    /* A is not present */
  {
    stats [STATUS] = ERROR_A_NOT_PRESENT ;
    _debug0 ("colamd: A not present\n") ;
    return (_false) ;
  }

  if (p == null)    /* p is not present */
  {
    stats [STATUS] = ERROR_P_NOT_PRESENT ;
    _debug0 ("colamd: p not present\n") ;
    return (_false) ;
  }

  if (n_row < 0)  /* n_row must be >= 0 */
  {
    stats [STATUS] = ERROR_NROW_NEGATIVE ;
    stats [INFO1] = n_row ;
    _debug0 ("colamd: nrow negative $n_row\n") ;
    return (_false) ;
  }

  if (n_col < 0)  /* n_col must be >= 0 */
  {
    stats [STATUS] = ERROR_NCOL_NEGATIVE ;
    stats [INFO1] = n_col ;
    _debug0 ("colamd: ncol negative $n_col\n") ;
    return (_false) ;
  }

  nnz = p [n_col] ;
  if (nnz < 0)  /* nnz must be >= 0 */
  {
    stats [STATUS] = ERROR_NNZ_NEGATIVE ;
    stats [INFO1] = nnz ;
    _debug0 ("colamd: number of entries negative $nnz\n") ;
    return (_false) ;
  }

  if (p [0] != 0)
  {
    stats [STATUS] = ERROR_P0_NONZERO  ;
    stats [INFO1] = p [0] ;
    _debug0 ("colamd: p[0] not zero ${p [0]}\n") ;
    return (_false) ;
  }

  /* === If no knobs, set default knobs =============================== */

  if (knobs == null)
  {
    setDefaults (default_knobs) ;
    knobs = default_knobs ;
  }

  aggressive = (knobs [AGGRESSIVE] != _false) ? 1 : 0;

  /* === Allocate the Row and Col arrays from array A ================= */

  ok = new Int32List.fromList([ _true ]) ;
  Col_size = _c (n_col, ok) ;      /* size of Col array of structs */
  Row_size = _r (n_row, ok) ;      /* size of Row array of structs */

  /* need = 2*nnz + n_col + Col_size + Row_size ; */
  need = _mult (nnz, 2, ok) ;
  need = _add (need, n_col, ok) ;
//    need = t_add (need, Col_size, ok) ;
//    need = t_add (need, Row_size, ok) ;

  if ((ok[0] == 0) || need > Alen.toInt() /*|| need > Int_MAX*/)
  {
    /* not enough space in array A to perform the ordering */
    stats [STATUS] = ERROR_A_TOO_SMALL ;
    stats [INFO1] = need ;
    stats [INFO2] = Alen ;
    _debug0 ("colamd: Need Alen >= $need, given only Alen = $Alen\n");
    return (_false) ;
  }

//    Alen -= Col_size + Row_size ;
  cols = new List<Col>(Col_size) ;  //A [Alen] ;
  rows = new List<Row>(Row_size) ;  //A [Alen + Col_size] ;

  for (int i = 0; i < Col_size; i++) {
    cols [i] = new Col() ;
  }
  for (int i = 0; i < Row_size; i++) {
    rows [i] = new Row() ;
  }

  /* === Construct the row and column data structures ================= */

  if (_initRowsCols (n_row, n_col, rows, cols, A, p, stats) == 0)
  {
    /* input matrix is invalid */
    _debug0 ("colamd: Matrix invalid\n") ;
    return (_false) ;
  }

  /* === Initialize scores, kill dense rows/columns =================== */

  _initScoring (n_row, n_col, rows, cols, A, p, knobs,
      n_row2, n_col2, max_deg) ;

  /* === Order the supercolumns ======================================= */

  ngarbage = _findOrdering (n_row, n_col, Alen, rows, cols, A, p,
      n_col2[0], max_deg[0], 2*nnz, aggressive) ;

  /* === Order the non-principal columns ============================== */

  _orderChildren (n_col, cols, p) ;

  /* === Return statistics in stats =================================== */

  stats [DENSE_ROW] = n_row - n_row2[0] ;
  stats [DENSE_COL] = n_col - n_col2[0] ;
  stats [DEFRAG_COUNT] = ngarbage ;
  _debug0 ("colamd: done.\n") ;
  return (_true) ;
}


/* ========================================================================== */
/* === colamd_report ======================================================== */
/* ========================================================================== */

void colamdReport(List<int> stats)
{
  _printReport ("colamd", stats) ;
}


/* ========================================================================== */
/* === symamd_report ======================================================== */
/* ========================================================================== */

void symamdReport(List<int> stats)
{
  _printReport ("symamd", stats) ;
}



/* ========================================================================== */
/* === NON-USER-CALLABLE ROUTINES: ========================================== */
/* ========================================================================== */

/* There are no user-callable routines beyond this point in the file */


/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

/// Takes the column form of the matrix in A and creates the row form of the
/// matrix.  Also, row and column attributes are stored in the Col and Row
/// structs.  If the columns are un-sorted or contain duplicate row indices,
/// this routine will also sort and remove duplicate row indices from the
/// column form of the matrix.  Returns FALSE if the matrix is invalid,
/// TRUE otherwise.  Not user-callable.
///
/// @param n_row number of rows of A
/// @param n_col number of columns of A
/// @param Row of size n_row+1
/// @param Col of size n_col+1
/// @param A row indices of A, of size Alen
/// @param p pointers to columns in A, of size n_col+1
/// @param stats colamd statistics
/// @return TRUE if OK, or FALSE otherwise
int _initRowsCols(int n_row, int n_col, List<Row> rows,
    List<Col> cols, List<int> A, List<int> p, List<int> stats)
{
  /* === Local variables ============================================== */

  int col ;    // a column index
  int row ;    // a row index
  int cp ;    // a column pointer
  int cp_end ;    // a pointer to the end of a column
  int rp ;    // a row pointer
  int rp_end ;    // a pointer to the end of a row
  int last_row ;    // previous row

  /* === Initialize columns, and check column pointers ================ */

  for (col = 0 ; col < n_col ; col++)
  {
    cols [col].start = p [col] ;
    cols [col].length = p [col+1] - p [col] ;

    if (cols [col].length < 0)
    {
      /* column pointers must be non-decreasing */
      stats [STATUS] = ERROR_COL_LENGTH_NEGATIVE ;
      stats [INFO1] = col ;
      stats [INFO2] = cols [col].length ;
      _debug0 ("colamd: col $col length ${cols [col].length} < 0\n") ;
      return (_false) ;
    }

    cols [col].thickness = 1 ;
    cols [col].score = 0 ;
    cols [col].prev = _empty ;
    cols [col].degree_next = _empty ;
  }

  /* p [0..n_col] no longer needed, used as "head" in subsequent routines */

  /* === Scan columns, compute row degrees, and check row indices ===== */

  stats [INFO3] = 0 ;  /* number of duplicate or unsorted row indices*/

  for (row = 0 ; row < n_row ; row++)
  {
    rows [row].length = 0 ;
    rows [row].mark = -1 ;
  }

  for (col = 0 ; col < n_col ; col++)
  {
    last_row = -1 ;

    cp = p [col] ;
    cp_end = p [col+1] ;

    while (cp < cp_end)
    {
      row = A [cp++] ;

      /* make sure row indices within range */
      if (row < 0 || row >= n_row)
      {
        stats [STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS ;
        stats [INFO1] = col ;
        stats [INFO2] = row ;
        stats [INFO3] = n_row ;
        _debug0 ("colamd: row $row col $col out of bounds\n") ;
        return (_false) ;
      }

      if (row <= last_row || rows [row].mark == col)
      {
        /* row index are unsorted or repeated (or both), thus col */
        /* is jumbled.  This is a notice, not an error condition. */
        stats [STATUS] = OK_BUT_JUMBLED ;
        stats [INFO1] = col ;
        stats [INFO2] = row ;
        stats [INFO3] ++ ;
        _debug1 ("colamd: row $row col $col unsorted/duplicate\n");
      }

      if (rows [row].mark != col)
      {
        rows [row].length++ ;
      }
      else
      {
        /* this is a repeated entry in the column, */
        /* it will be removed */
        cols [col].length-- ;
      }

      /* mark the row as having been seen in this column */
      rows [row].mark = col ;

      last_row = row ;
    }
  }

  /* === Compute row pointers ========================================= */

  /* row form of the matrix starts directly after the column */
  /* form of matrix in A */
  rows [0].start = p [n_col] ;
  rows [0].p = rows [0].start ;
  rows [0].mark = -1 ;
  for (row = 1 ; row < n_row ; row++)
  {
    rows [row].start = rows [row-1].start + rows [row-1].length ;
    rows [row].p = rows [row].start ;
    rows [row].mark = -1 ;
  }

  /* === Create row form ============================================== */

  if (stats [STATUS] == OK_BUT_JUMBLED)
  {
    /* if cols jumbled, watch for repeated row indices */
    for (col = 0 ; col < n_col ; col++)
    {
      cp = p [col] ;
      cp_end = p [col+1] ;

      while (cp < cp_end)
      {
        row = A [cp++] ;

        if (rows [row].mark != col)
        {
          A [rows [row].p] = col ;
          rows [row].p = rows [row].p + 1 ;
          rows [row].mark = col ;
        }
      }
    }
  }
  else
  {
    /* if cols not jumbled, we don't need the mark (this is faster) */
    for (col = 0 ; col < n_col ; col++)
    {
      cp = p [col] ;
      cp_end = p [col+1] ;
      while (cp < cp_end)
      {
        A [rows [A [cp]].p] = col ;
        rows [A [cp]].p = rows [A [cp]].p + 1 ;
        cp++ ;
      }
    }
  }

  /* === Clear the row marks and set row degrees ====================== */

  for (row = 0 ; row < n_row ; row++)
  {
    rows [row].mark = 0 ;
    rows [row].degree = rows [row].length ;
  }

  /* === See if we need to re-create columns ========================== */

  if (stats [STATUS] == OK_BUT_JUMBLED)
  {
    _debug0 ("colamd: reconstructing column form, matrix jumbled\n") ;

    if (!_ndebug)
    {
      /* make sure column lengths are correct */
      for (col = 0 ; col < n_col ; col++)
      {
        p [col] = cols [col].length ;
      }
      for (row = 0 ; row < n_row ; row++)
      {
        rp = rows [row].start ;
        rp_end = rp + rows [row].length ;
        while (rp < rp_end)
        {
          p [A [rp++]]-- ;
        }
      }
      for (col = 0 ; col < n_col ; col++)
      {
        _assert (p [col] == 0) ;
      }
      /* now p is all zero (different than when debugging is turned off) */
    } /* NDEBUG */

    /* === Compute col pointers ========================================= */

    /* col form of the matrix starts at A [0]. */
    /* Note, we may have a gap between the col form and the row */
    /* form if there were duplicate entries, if so, it will be */
    /* removed upon the first garbage collection */
    cols [0].start = 0 ;
    p [0] = cols [0].start ;
    for (col = 1 ; col < n_col ; col++)
    {
      /* note that the lengths here are for pruned columns, i.e. */
      /* no duplicate row indices will exist for these columns */
      cols [col].start = cols [col-1].start + cols [col-1].length ;
      p [col] = cols [col].start ;
    }

    /* === Re-create col form =========================================== */

    for (row = 0 ; row < n_row ; row++)
    {
      rp = rows [row].start ;
      rp_end = rp + rows [row].length ;
      while (rp < rp_end)
      {
        A [p [A [rp++]]++] = row ;
      }
    }
  }

  /* === Done.  Matrix is not (or no longer) jumbled ================== */

  return _true;
}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

/**
 * Kills dense or empty columns and rows, calculates an initial score for
 * each column, and places all columns in the degree lists.init_rows_cols
 *
 * @param n_row number of rows of A
 * @param n_col number of columns of A
 * @param Row of size n_row+1
 * @param Col of size n_col+1
 * @param A column form and row form of A
 * @param head of size n_col+1
 * @param knobs parameters
 * @param p_n_row2 size 1, number of non-dense, non-empty rows
 * @param p_n_col2 size 1, number of non-dense, non-empty columns
 * @param p_max_deg size 1, maximum row degree
 */
void _initScoring (int n_row, int n_col, List<Row> rows,
    List<Col> cols, Int32List A, Int32List head, Float64List knobs,
    Int32List p_n_row2, Int32List p_n_col2, Int32List p_max_deg)
{
  /* === Local variables ============================================== */

  int c ;                // a column index
  int r, row ;           // a row index
  int cp ;               // a column pointer
  int deg ;              // degree of a row or column
  int cp_end ;           // a pointer to the end of a column
  int new_cp ;           // new column pointer
  int col_length ;       // length of pruned column
  int score ;            // current column score
  int n_col2 ;           // number of non-dense, non-empty columns
  int n_row2 ;           // number of non-dense, non-empty rows
  int dense_row_count ;  // remove rows with more entries than this
  int dense_col_count ;  // remove cols with more entries than this
  int min_score ;        // smallest column score
  int max_deg ;          // maximum row degree
  int next_col ;         // Used to add to degree list.

  int debug_count = 0 ;  // debug only.

  /* === Extract knobs ================================================ */

  /* Note: if knobs contains a NaN, this is undefined: */
  if (knobs [DENSE_ROW] < 0)
  {
    /* only remove completely dense rows */
    dense_row_count = n_col-1 ;
  }
  else
  {
    dense_row_count = _denseDegree (knobs [DENSE_ROW], n_col) ;
  }
  if (knobs [DENSE_COL] < 0)
  {
    /* only remove completely dense columns */
    dense_col_count = n_row-1 ;
  }
  else
  {
    dense_col_count = _denseDegree (knobs [DENSE_COL],
        _min (n_row, n_col)) ;
  }

  _debug1 ("colamd: densecount: $dense_row_count $dense_col_count\n") ;
  max_deg = 0 ;
  n_col2 = n_col ;
  n_row2 = n_row ;

  /* === Kill empty columns =========================================== */

  /* Put the empty columns at the end in their natural order, so that LU */
  /* factorization can proceed as far as possible. */
  for (c = n_col-1 ; c >= 0 ; c--)
  {
    deg = cols [c].length ;
    if (deg == 0)
    {
      /* this is a empty column, kill and order it last */
      cols [c].order = --n_col2 ;
      _killPrincipalCol (cols, c) ;
    }
  }
  _debug1 ("colamd: null columns killed: ${n_col - n_col2}\n") ;

  /* === Kill dense columns =========================================== */

  /* Put the dense columns at the end, in their natural order */
  for (c = n_col-1 ; c >= 0 ; c--)
  {
    /* skip any dead columns */
    if (_colIsDead (cols, c))
    {
      continue ;
    }
    deg = cols [c].length ;
    if (deg > dense_col_count)
    {
      /* this is a dense column, kill and order it last */
      cols [c].order = --n_col2 ;
      /* decrement the row degrees */
      cp = cols [c].start ;
      cp_end = cp + cols [c].length ;
      while (cp < cp_end)
      {
        rows [A [cp]].degree = rows [A [cp]].degree - 1 ;
        cp++ ;
      }
      _killPrincipalCol (cols, c) ;
    }
  }
  _debug1 ("colamd: Dense and null columns killed: ${n_col - n_col2}\n") ;

  /* === Kill dense and empty rows ==================================== */

  for (r = 0 ; r < n_row ; r++)
  {
    deg = rows [r].degree ;
    _assert (deg >= 0 && deg <= n_col) ;
    if (deg > dense_row_count || deg == 0)
    {
      /* kill a dense or empty row */
      _killRow (rows, r) ;
      --n_row2 ;
    }
    else
    {
      /* keep track of max degree of remaining rows */
      max_deg = _max (max_deg, deg) ;
    }
  }
  _debug1 ("colamd: Dense and null rows killed: ${n_row - n_row2}\n") ;

  /* === Compute initial column scores ================================ */

  /* At this point the row degrees are accurate.  They reflect the number */
  /* of "live" (non-dense) columns in each row.  No empty rows exist. */
  /* Some "live" columns may contain only dead rows, however.  These are */
  /* pruned in the code below. */

  /* now find the initial matlab score for each column */
  for (c = n_col-1 ; c >= 0 ; c--)
  {
    /* skip dead column */
    if (_colIsDead (cols, c))
    {
      continue ;
    }
    score = 0 ;
    cp = cols [c].start ;
    new_cp = cp ;
    cp_end = cp + cols [c].length ;
    while (cp < cp_end)
    {
      /* get a row */
      row = A [cp++] ;
      /* skip if dead */
      if (_rowIsDead (rows, row))
      {
        continue ;
      }
      /* compact the column */
      A [new_cp++] = row ;
      /* add row's external degree */
      score += rows [row].degree - 1 ;
      /* guard against integer overflow */
      score = _min (score, n_col) ;
    }
    /* determine pruned column length */
    col_length = (new_cp - cols [c].start) ;
    if (col_length == 0)
    {
      /* a newly-made null column (all rows in this col are "dense" */
      /* and have already been killed) */
      _debug2 ("Newly null killed: $c\n") ;
      cols [c].order = --n_col2 ;
      _killPrincipalCol (cols, c) ;
    }
    else
    {
      /* set column length and set score */
      _assert (score >= 0) ;
      _assert (score <= n_col) ;
      cols [c].length = col_length ;
      cols [c].score = score ;
    }
  }
  _debug1 ("colamd: Dense, null, and newly-null columns killed: ${n_col-n_col2}\n") ;

  /* At this point, all empty rows and columns are dead.  All live columns */
  /* are "clean" (containing no dead rows) and simplicial (no supercolumns */
  /* yet).  Rows may contain dead columns, but all live rows contain at */
  /* least one live column. */

  if (!_ndebug)
  {
    _debugStructures (n_row, n_col, rows, cols, A, n_col2) ;
  }

  /* === Initialize degree lists ========================================== */

  if (!_ndebug)
  {
    debug_count = 0 ;
  }

  /* clear the hash buckets */
  for (c = 0 ; c <= n_col ; c++)
  {
    head [c] = _empty ;
  }
  min_score = n_col ;
  /* place in reverse order, so low column indices are at the front */
  /* of the lists.  This is to encourage natural tie-breaking */
  for (c = n_col-1 ; c >= 0 ; c--)
  {
    /* only add principal columns to degree lists */
    if (_colIsAlive (cols, c))
    {
      _debug4 ("place $c score ${cols [c].score} minscore $min_score ncol $n_col\n") ;

      /* === Add columns score to DList =============================== */

      score = cols [c].score ;

      _assert (min_score >= 0) ;
      _assert (min_score <= n_col) ;
      _assert (score >= 0) ;
      _assert (score <= n_col) ;
      _assert (head [score] >= _empty) ;

      /* now add this column to dList at proper score location */
      next_col = head [score] ;
      cols [c].prev = _empty ;
      cols [c].degree_next = next_col ;

      /* if there already was a column with the same score, set its */
      /* previous pointer to this new column */
      if (next_col != _empty)
      {
        cols [next_col].prev = c ;
      }
      head [score] = c ;

      /* see if this score is less than current min */
      min_score = _min (min_score, score) ;

      if (!_ndebug)
      {
        debug_count++ ;
      }

    }
  }

  if (!_ndebug)
  {
    _debug1 ("colamd: Live cols $debug_count out of $n_col, non-princ: ${n_col-debug_count}\n") ;
    _assert (debug_count == n_col2) ;
    _debugDegLists (n_row, n_col, rows, cols, head, min_score, n_col2, max_deg) ;
  } /* NDEBUG */

  /* === Return number of remaining columns, and max row degree ======= */

  p_n_col2[0] = n_col2 ;
  p_n_row2[0] = n_row2 ;
  p_max_deg[0] = max_deg ;
}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

///
/// Order the principal columns of the supercolumn form of the matrix
/// (no supercolumns on input).  Uses a minimum approximate column minimum
/// degree ordering method.  Not user-callable.
///
/// @param n_row number of rows of A
/// @param n_col number of columns of A
/// @param Alen size of A, 2*nnz + n_col or larger
/// @param Row of size n_row+1
/// @param Col of size n_col+1
/// @param A column form and row form of A
/// @param head of size n_col+1
/// @param n_col2 Remaining columns to order
/// @param max_deg Maximum row degree
/// @param pfree index of first free slot (2*nnz on entry)
/// @param aggressive
/// @return the number of garbage collections
int _findOrdering (int n_row, int n_col, int Alen,
    List<Row> rows, List<Col> cols, List<int> A, List<int> head, int n_col2,
    int max_deg, int pfree, int aggressive)
{
  /* === Local variables ============================================== */

  int k ;                 // current pivot ordering step
  int pivot_col ;         // current pivot column
  int cp ;                // a column pointer
  int rp ;                // a row pointer
  int pivot_row ;         // current pivot row
  int new_cp ;            // modified column pointer
  int new_rp ;            // modified row pointer
  int pivot_row_start ;   // pointer to start of pivot row
  int pivot_row_degree ;  // number of columns in pivot row
  int pivot_row_length ;  // number of supercolumns in pivot row
  int pivot_col_score ;   // score of pivot column
  int needed_memory ;     // free space needed for pivot row
  int cp_end ;            // pointer to the end of a column
  int rp_end ;            // pointer to the end of a row
  int row ;               // a row index
  int col ;               // a column index
  int max_score ;         // maximum possible score
  int cur_score ;         // score of current column
  /* FIXME unsigned */ int hash ; // hash value for supernode detection
  int head_column ;       // head of hash bucket
  int first_col ;         // first column in hash bucket
  int tag_mark ;          // marker value for mark array
  int row_mark ;          // Row [row].shared2.mark
  int set_difference ;    // set difference size of row with pivot row
  int min_score ;         // smallest column score
  int col_thickness ;     // "thickness" (no. of columns in a supercol)
  int max_mark ;          // maximum value of tag_mark
  int pivot_col_thickness ; // number of columns represented by pivot col
  int prev_col ;          // Used by Dlist operations.
  int next_col ;          // Used by Dlist operations.
  int ngarbage ;          // number of garbage collections performed

  int debug_d ;           // debug loop counter
  int debug_step = 0 ;    // debug loop counter

  /* === Initialization and clear mark ================================ */

  max_mark = double.MAX_FINITE.toInt();//Int_MAX - n_col ;
  tag_mark = _clearMark (0, max_mark, n_row, rows) ;
  min_score = 0 ;
  ngarbage = 0 ;
  _debug1 ("colamd: Ordering, n_col2=$n_col2\n") ;

  /* === Order the columns ============================================ */

  for (k = 0 ; k < n_col2 ; /* 'k' is incremented below */)
  {

    if (!_ndebug)
    {
      if (debug_step % 100 == 0)
      {
        _debug2 ("\n...       Step k: $k out of n_col2: $n_col2\n") ;
      }
      else
      {
        _debug3 ("\n----------Step k: $k out of n_col2: $n_col2\n") ;
      }
      debug_step++ ;
      _debugDegLists (n_row, n_col, rows, cols, head,
          min_score, n_col2-k, max_deg) ;
      _debugMatrix (n_row, n_col, rows, cols, A) ;
    } /* NDEBUG */

    /* === Select pivot column, and order it ============================ */

    /* make sure degree list isn't empty */
    _assert (min_score >= 0) ;
    _assert (min_score <= n_col) ;
    _assert (head [min_score] >= _empty) ;

    if (!_ndebug)
    {
      for (debug_d = 0 ; debug_d < min_score ; debug_d++)
      {
        _assert (head [debug_d] == _empty) ;
      }
    } /* NDEBUG */

    /* get pivot column from head of minimum degree list */
    while (head [min_score] == _empty && min_score < n_col)
    {
      min_score++ ;
    }
    pivot_col = head [min_score] ;
    _assert (pivot_col >= 0 && pivot_col <= n_col) ;
    next_col = cols [pivot_col].degree_next ;
    head [min_score] = next_col ;
    if (next_col != _empty)
    {
      cols [next_col].prev = _empty ;
    }

    _assert (_colIsAlive (cols, pivot_col)) ;

    /* remember score for defrag check */
    pivot_col_score = cols [pivot_col].score ;

    /* the pivot column is the kth column in the pivot order */
    cols [pivot_col].order = k ;

    /* increment order count by column thickness */
    pivot_col_thickness = cols [pivot_col].thickness ;
    k += pivot_col_thickness ;
    _assert (pivot_col_thickness > 0) ;
    _debug3 ("Pivot col: $pivot_col thick $pivot_col_thickness\n") ;

    /* === Garbage_collection, if necessary ============================= */

    needed_memory = _min (pivot_col_score, n_col - k) ;
    if (pfree + needed_memory >= Alen)
    {
      pfree = _garbageCollection (n_row, n_col, rows, cols, A, pfree) ;
      ngarbage++ ;
      /* after garbage collection we will have enough */
      _assert (pfree + needed_memory < Alen) ;
      /* garbage collection has wiped out the Row[].shared2.mark array */
      tag_mark = _clearMark (0, max_mark, n_row, rows) ;

      if (!_ndebug)
      {
        _debugMatrix (n_row, n_col, rows, cols, A) ;
      } /* NDEBUG */
    }

    /* === Compute pivot row pattern ==================================== */

    /* get starting location for this new merged row */
    pivot_row_start = pfree ;

    /* initialize new row counts to zero */
    pivot_row_degree = 0 ;

    /* tag pivot column as having been visited so it isn't included */
    /* in merged pivot row */
    cols [pivot_col].thickness = -pivot_col_thickness ;

    /* pivot row is the union of all rows in the pivot column pattern */
    cp = cols [pivot_col].start ;
    cp_end = cp + cols [pivot_col].length ;
    while (cp < cp_end)
    {
      /* get a row */
      row = A [cp++] ;
      _debug4 ("Pivot col pattern ${_rowIsAlive (rows, row) ? 1 : 0} $row\n") ;
      /* skip if row is dead */
      if (_rowIsAlive (rows, row))
      {
        rp = rows [row].start ;
        rp_end = rp + rows [row].length ;
        while (rp < rp_end)
        {
          /* get a column */
          col = A [rp++] ;
          /* add the column, if alive and untagged */
          col_thickness = cols [col].thickness ;
          if (col_thickness > 0 && _colIsAlive (cols, col))
          {
            /* tag column in pivot row */
            cols [col].thickness = -col_thickness ;
            _assert (pfree < Alen) ;
            /* place column in pivot row */
            A [pfree++] = col ;
            pivot_row_degree += col_thickness ;
          }
        }
      }
    }

    /* clear tag on pivot column */
    cols [pivot_col].thickness = pivot_col_thickness ;
    max_deg = _max (max_deg, pivot_row_degree) ;

    if (!_ndebug)
    {
      _debug3 ("check2\n") ;
      _debugMark (n_row, rows, tag_mark, max_mark) ;
    } /* NDEBUG */

    /* === Kill all rows used to construct pivot row ==================== */

    /* also kill pivot row, temporarily */
    cp = cols [pivot_col].start ;
    cp_end = cp + cols [pivot_col].length ;
    while (cp < cp_end)
    {
      /* may be killing an already dead row */
      row = A [cp++] ;
      _debug3 ("Kill row in pivot col: $row\n") ;
      _killRow (rows, row) ;
    }

    /* === Select a row index to use as the new pivot row =============== */

    pivot_row_length = pfree - pivot_row_start ;
    if (pivot_row_length > 0)
    {
      /* pick the "pivot" row arbitrarily (first row in col) */
      pivot_row = A [cols [pivot_col].start] ;
      _debug3 ("Pivotal row is $pivot_row\n") ;
    }
    else
    {
      /* there is no pivot row, since it is of zero length */
      pivot_row = _empty ;
      _assert (pivot_row_length == 0) ;
    }
    _assert (cols [pivot_col].length > 0 || pivot_row_length == 0) ;

    /* === Approximate degree computation =============================== */

    /* Here begins the computation of the approximate degree.  The column */
    /* score is the sum of the pivot row "length", plus the size of the */
    /* set differences of each row in the column minus the pattern of the */
    /* pivot row itself.  The column ("thickness") itself is also */
    /* excluded from the column score (we thus use an approximate */
    /* external degree). */

    /* The time taken by the following code (compute set differences, and */
    /* add them up) is proportional to the size of the data structure */
    /* being scanned - that is, the sum of the sizes of each column in */
    /* the pivot row.  Thus, the amortized time to compute a column score */
    /* is proportional to the size of that column (where size, in this */
    /* context, is the column "length", or the number of row indices */
    /* in that column).  The number of row indices in a column is */
    /* monotonically non-decreasing, from the length of the original */
    /* column on input to colamd. */

    /* === Compute set differences ====================================== */

    _debug3 ("** Computing set differences phase. **\n") ;

    /* pivot row is currently dead - it will be revived later. */

    _debug3 ("Pivot row: ") ;
    /* for each column in pivot row */
    rp = pivot_row_start ;
    rp_end = rp + pivot_row_length ;
    while (rp < rp_end)
    {
      col = A [rp++] ;
      _assert (_colIsAlive (cols, col) && col != pivot_col) ;
      _debug3 ("Col: $col\n") ;

      /* clear tags used to construct pivot row pattern */
      col_thickness = -cols [col].thickness ;
      _assert (col_thickness > 0) ;
      cols [col].thickness = col_thickness ;

      /* === Remove column from degree list =========================== */

      cur_score = cols [col].score ;
      prev_col = cols [col].prev ;
      next_col = cols [col].degree_next ;
      _assert (cur_score >= 0) ;
      _assert (cur_score <= n_col) ;
      _assert (cur_score >= _empty) ;
      if (prev_col == _empty)
      {
        head [cur_score] = next_col ;
      }
      else
      {
        cols [prev_col].degree_next = next_col ;
      }
      if (next_col != _empty)
      {
        cols [next_col].prev = prev_col ;
      }

      /* === Scan the column ========================================== */

      cp = cols [col].start ;
      cp_end = cp + cols [col].length ;
      while (cp < cp_end)
      {
        /* get a row */
        row = A [cp++] ;
        row_mark = rows [row].mark ;
        /* skip if dead */
        if (_rowIsMarkedDead (row_mark))
        {
          continue ;
        }
        _assert (row != pivot_row) ;
        set_difference = row_mark - tag_mark ;
        /* check if the row has been seen yet */
        if (set_difference < 0)
        {
          _assert (rows [row].degree <= max_deg) ;
          set_difference = rows [row].degree ;
        }
        /* subtract column thickness from this row's set difference */
        set_difference -= col_thickness ;
        _assert (set_difference >= 0) ;
        /* absorb this row if the set difference becomes zero */
        if (set_difference == 0 && aggressive != 0)
        {
          _debug3 ("aggressive absorption. Row: $row\n") ;
          _killRow (rows, row) ;
        }
        else
        {
          /* save the new mark */
          rows [row].mark = set_difference + tag_mark ;
        }
      }
    }

    if (!_ndebug)
    {
      _debugDegLists (n_row, n_col, rows, cols, head,
          min_score, n_col2-k-pivot_row_degree, max_deg) ;
    } /* NDEBUG */

    /* === Add up set differences for each column ======================= */

    _debug3 ("** Adding set differences phase. **\n") ;

    /* for each column in pivot row */
    rp = pivot_row_start ;
    rp_end = rp + pivot_row_length ;
    while (rp < rp_end)
    {
      /* get a column */
      col = A [rp++] ;
      _assert (_colIsAlive (cols, col) && col != pivot_col) ;
      hash = 0 ;
      cur_score = 0 ;
      cp = cols [col].start ;
      /* compact the column */
      new_cp = cp ;
      cp_end = cp + cols [col].length ;

      _debug4 ("Adding set diffs for Col: $col.\n") ;

      while (cp < cp_end)
      {
        /* get a row */
        row = A [cp++] ;
        _assert(row >= 0 && row < n_row) ;
        row_mark = rows [row].mark ;
        /* skip if dead */
        if (_rowIsMarkedDead (row_mark))
        {
          _debug4 (" Row $row, dead\n") ;
          continue ;
        }
        _debug4 (" Row $row, set diff ${row_mark-tag_mark}\n");
        _assert (row_mark >= tag_mark) ;
        /* compact the column */
        A [new_cp++] = row ;
        /* compute hash function */
        hash += row ;
        /* add set difference */
        cur_score += row_mark - tag_mark ;
        /* integer overflow... */
        cur_score = _min (cur_score, n_col) ;
      }

      /* recompute the column's length */
      cols [col].length = (new_cp - cols [col].start) ;

      /* === Further mass elimination ================================= */

      if (cols [col].length == 0)
      {
        _debug4 ("further mass elimination. Col: $col\n") ;
        /* nothing left but the pivot row in this column */
        _killPrincipalCol (cols, col) ;
        pivot_row_degree -= cols [col].thickness ;
        _assert (pivot_row_degree >= 0) ;
        /* order it */
        cols [col].order = k ;
        /* increment order count by column thickness */
        k += cols [col].thickness ;
      }
      else
      {
        /* === Prepare for supercolumn detection ==================== */

        _debug4 ("Preparing supercol detection for Col: $col.\n") ;

        /* save score so far */
        cols [col].score = cur_score ;

        /* add column to hash table, for supercolumn detection */
        hash %= n_col + 1 ;

        _debug4 (" Hash = $hash, n_col = $n_col.\n") ;
        _assert (hash.toInt() <= n_col) ;

        head_column = head [hash] ;
        if (head_column > _empty)
        {
          /* degree list "hash" is non-empty, use prev (shared3) of */
          /* first column in degree list as head of hash bucket */
          first_col = cols [head_column].headhash ;
          cols [head_column].headhash = col ;
        }
        else
        {
          /* degree list "hash" is empty, use head as hash bucket */
          first_col = - (head_column + 2) ;
          head [hash] = - (col + 2) ;
        }
        cols [col].hash_next = first_col ;

        /* save hash function in Col [col].shared3.hash */
        cols [col].hash = hash.toInt() ;
        _assert (_colIsAlive (cols, col)) ;
      }
    }

    /* The approximate external column degree is now computed.  */

    /* === Supercolumn detection ======================================== */

    _debug3 ("** Supercolumn detection phase. **\n") ;

    if (!_ndebug)
    {
      _detectSuperCols (n_col, rows,
          cols, A, head, pivot_row_start, pivot_row_length) ;
    } else {
      _detectSuperCols (0, null, cols, A, head, pivot_row_start, pivot_row_length) ;
    }

    /* === Kill the pivotal column ====================================== */

    _killPrincipalCol (cols, pivot_col) ;

    /* === Clear mark =================================================== */

    tag_mark = _clearMark (tag_mark+max_deg+1, max_mark, n_row, rows) ;

    if (!_ndebug)
    {
      _debug3 ("check3\n") ;
      _debugMark (n_row, rows, tag_mark, max_mark) ;
    } /* NDEBUG */

    /* === Finalize the new pivot row, and column scores ================ */

    _debug3 ("** Finalize scores phase. **\n") ;

    /* for each column in pivot row */
    rp = pivot_row_start ;
    /* compact the pivot row */
    new_rp = rp ;
    rp_end = rp + pivot_row_length ;
    while (rp < rp_end)
    {
      col = A [rp++] ;
      /* skip dead columns */
      if (_colIsDead (cols, col))
      {
        continue ;
      }
      A [new_rp++] = col ;
      /* add new pivot row to column */
      A [cols [col].start + (cols [col].length++)] = pivot_row ;

      /* retrieve score so far and add on pivot row's degree. */
      /* (we wait until here for this in case the pivot */
      /* row's degree was reduced due to mass elimination). */
      cur_score = cols [col].score + pivot_row_degree ;

      /* calculate the max possible score as the number of */
      /* external columns minus the 'k' value minus the */
      /* columns thickness */
      max_score = n_col - k - cols [col].thickness ;

      /* make the score the external degree of the union-of-rows */
      cur_score -= cols [col].thickness ;

      /* make sure score is less or equal than the max score */
      cur_score = _min (cur_score, max_score) ;
      _assert (cur_score >= 0) ;

      /* store updated score */
      cols [col].score = cur_score ;

      /* === Place column back in degree list ========================= */

      _assert (min_score >= 0) ;
      _assert (min_score <= n_col) ;
      _assert (cur_score >= 0) ;
      _assert (cur_score <= n_col) ;
      _assert (head [cur_score] >= _empty) ;
      next_col = head [cur_score] ;
      cols [col].degree_next = next_col ;
      cols [col].prev = _empty ;
      if (next_col != _empty)
      {
        cols [next_col].prev = col ;
      }
      head [cur_score] = col ;

      /* see if this score is less than current min */
      min_score = _min (min_score, cur_score) ;

    }

    if (!_ndebug)
    {
      _debugDegLists (n_row, n_col, rows, cols, head,
          min_score, n_col2-k, max_deg) ;
    } /* NDEBUG */

    /* === Resurrect the new pivot row ================================== */

    if (pivot_row_degree > 0)
    {
      /* update pivot row length to reflect any cols that were killed */
      /* during super-col detection and mass elimination */
      rows [pivot_row].start  = pivot_row_start ;
      rows [pivot_row].length = new_rp - pivot_row_start ;
      _assert (rows [pivot_row].length > 0) ;
      rows [pivot_row].degree = pivot_row_degree ;
      rows [pivot_row].mark = 0 ;
      /* pivot row is no longer dead */

      _debug1 ("Resurrect Pivot_row $pivot_row deg: $pivot_row_degree\n") ;
    }
  }

  /* === All principal columns have now been ordered ================== */

  return (ngarbage) ;
}


/* ========================================================================== */
/* === order_children ======================================================= */
/* ========================================================================== */

/// The find_ordering routine has ordered all of the principal columns (the
/// representatives of the supercolumns).  The non-principal columns have not
/// yet been ordered.  This routine orders those columns by walking up the
/// parent tree (a column is a child of the column which absorbed it).  The
/// final permutation vector is then placed in p [0 ... n_col-1], with p [0]
/// being the first column, and p [n_col-1] being the last.  It doesn't look
/// like it at first glance, but be assured that this routine takes time linear
/// in the number of columns.  Although not immediately obvious, the time
/// taken by this routine is O (n_col), that is, linear in the number of
/// columns.  Not user-callable.
///
/// @param n_col number of columns of A
/// @param Col of size n_col+1
/// @param p p [0 ... n_col-1] is the column permutation
void _orderChildren (int n_col, List<Col> cols, List<int> p)
{
  /* === Local variables ============================================== */

  int parent ;     /* index of column's parent */
  int order ;      /* column's order */

  /* === Order each non-principal column ============================== */

  for (int i = 0 ; i < n_col ; i++)
  {
    /* find an un-ordered non-principal column */
    _assert (_colIsDead (cols, i)) ;
    if (!_colIsDeadPrincipal (cols, i) && cols [i].order == _empty)
    {
      parent = i ;
      /* once found, find its principal parent */
      do
      {
        parent = cols [parent].parent ;
      } while (!_colIsDeadPrincipal (cols, parent)) ;

      /* now, order all un-ordered non-principal columns along path */
      /* to this parent.  collapse tree at the same time */
      int c = i ;
      /* get order of parent */
      order = cols [parent].order ;

      do
      {
        _assert (cols [c].order == _empty) ;

        /* order this column */
        cols [c].order = order++ ;
        /* collaps tree */
        cols [c].parent = parent ;

        /* get immediate parent of this column */
        c = cols [c].parent ;

        /* continue until we hit an ordered column.  There are */
        /* guarranteed not to be anymore unordered columns */
        /* above an ordered column */
      } while (cols [c].order == _empty) ;

      /* re-order the super_col parent to largest order for this group */
      cols [parent].order = order ;
    }
  }

  /* === Generate the permutation ===================================== */

  for (int c = 0 ; c < n_col ; c++)
  {
    p [cols [c].order] = c ;
  }
}


/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

/// Detects supercolumns by finding matches between columns in the hash buckets.
/// Check amongst columns in the set A [row_start ... row_start + row_length-1].
/// The columns under consideration are currently *not* in the degree lists,
/// and have already been placed in the hash buckets.
///
/// The hash bucket for columns whose hash function is equal to h is stored
/// as follows:
///
/// if head [h] is >= 0, then head [h] contains a degree list, so:
///
///     head [h] is the first column in degree bucket h.
///     Col [head [h]].headhash gives the first column in hash bucket h.
///
/// otherwise, the degree list is empty, and:
///
///     -(head [h] + 2) is the first column in hash bucket h.
///
/// For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
/// column" pointer.  Col [c].shared3.hash is used instead as the hash number
/// for that column.  The value of Col [c].shared4.hash_next is the next column
/// in the same hash bucket.
///
/// Assuming no, or "few" hash collisions, the time taken by this routine is
/// linear in the sum of the sizes (lengths) of each column whose score has
/// just been computed in the approximate degree computation.
/// Not user-callable.
///
/// @param n_col number of columns of A
/// @param Row of size n_row+1
/// @param Col of size n_col+1
/// @param A row indices of A
/// @param head head of degree lists and hash buckets
/// @param row_start pointer to set of columns to check
/// @param row_length number of columns to check
void _detectSuperCols (int n_col, List<Row> rows,
  List<Col> cols, List<int> A, List<int> head, int row_start, int row_length)
{
  /* === Local variables ============================================== */

  int hash ;         // hash value for a column
  int rp ;           // pointer to a row
  int c ;            // a column index
  int super_c ;      // column index of the column to absorb into
  int cp1 ;          // column pointer for column super_c
  int cp2 ;          // column pointer for column c
  int length ;       // length of column super_c
  int prev_c ;       // column preceding c in hash bucket
  int i ;            // loop counter
  int rp_end ;       // pointer to the end of the row
  int col ;          // a column index in the row to check
  int head_column ;  // first column in hash bucket or degree list
  int first_col ;    // first column in hash bucket

  /* === Consider each column in the row ============================== */

  rp = row_start ;
  rp_end = rp + row_length ;
  while (rp < rp_end)
  {
    col = A [rp++] ;
    if (_colIsDead (cols, col))
    {
      continue ;
    }

    /* get hash number for this column */
    hash = cols [col].hash ;
    _assert (hash <= n_col) ;

    /* === Get the first column in this hash bucket ===================== */

    head_column = head [hash] ;
    if (head_column > _empty)
    {
      first_col = cols [head_column].headhash ;
    }
    else
    {
      first_col = - (head_column + 2) ;
    }

    /* === Consider each column in the hash bucket ====================== */

    for (super_c = first_col ; super_c != _empty ;
        super_c = cols [super_c].hash_next)
    {
      _assert (_colIsAlive (cols, super_c)) ;
      _assert (cols [super_c].hash == hash) ;
      length = cols [super_c].length ;

      /* prev_c is the column preceding column c in the hash bucket */
      prev_c = super_c ;

      /* === Compare super_c with all columns after it ================ */

      for (c = cols [super_c].hash_next ;
          c != _empty ; c = cols [c].hash_next)
      {
        _assert (c != super_c) ;
        _assert (_colIsAlive (cols, c)) ;
        _assert (cols [c].hash == hash) ;

        /* not identical if lengths or scores are different */
        if (cols [c].length != length ||
            cols [c].score != cols [super_c].score)
        {
          prev_c = c ;
          continue ;
        }

        /* compare the two columns */
        cp1 = cols [super_c].start ;
        cp2 = cols [c].start ;

        for (i = 0 ; i < length ; i++)
        {
          /* the columns are "clean" (no dead rows) */
          _assert (_rowIsAlive (rows, A [cp1]))  ;
          _assert (_rowIsAlive (rows, A [cp2]))  ;
          /* row indices will same order for both supercols, */
          /* no gather scatter nessasary */
          if (A [cp1++] != A [cp2++])
          {
            break ;
          }
        }

        /* the two columns are different if the for-loop "broke" */
        if (i != length)
        {
          prev_c = c ;
          continue ;
        }

        /* === Got it!  two columns are identical =================== */

        _assert (cols [c].score == cols [super_c].score) ;

        cols [super_c].thickness = cols [super_c].thickness + cols [c].thickness ;
        cols [c].parent = super_c ;
        _killNonPrincipalCol (cols, c) ;
        /* order c later, in order_children() */
        cols [c].order = _empty ;
        /* remove c from hash bucket */
        cols [prev_c].hash_next = cols [c].hash_next ;
      }
    }

    /* === Empty this hash bucket ======================================= */

    if (head_column > _empty)
    {
      /* corresponding degree list "hash" is not empty */
      cols [head_column].headhash = _empty ;
    }
    else
    {
      /* corresponding degree list "hash" is empty */
      head [hash] = _empty ;
    }
  }
}


/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

/// Defragments and compacts columns and rows in the workspace A.  Used when
/// all avaliable memory has been used while performing row merging.  Returns
/// the index of the first free position in A, after garbage collection.  The
/// time taken by this routine is linear in the size of the array A, which is
/// itself linear in the number of nonzeros in the input matrix.
/// Not user-callable.
///
/// @param n_row number of rows
/// @param n_col number of columns
/// @param Row row info
/// @param Col column info
/// @param A A [0 ... Alen-1] holds the matrix
/// @param pfree
/// @return the new value of pfree
int _garbageCollection (int n_row, int n_col,
    List<Row> rows, List<Col> cols, List<int> A, int pfree)
{
  /* === Local variables ============================================== */

  int psrc ;    /* source pointer */
  int pdest ;   /* destination pointer */
  int j ;       /* counter */
  int r ;       /* a row index */
  int c ;       /* a column index */
  int length ;  /* length of a row or column */

  int debug_rows = 0 ;
  if (!_ndebug)
  {
    _debug2 ("Defrag..\n") ;
    for (psrc = 0 ; psrc < pfree ; psrc++) _assert (A [psrc] >= 0) ;
    debug_rows = 0 ;
  }

  /* === Defragment the columns ======================================= */

  pdest = 0 ;
  for (c = 0 ; c < n_col ; c++)
  {
    if (_colIsAlive (cols, c))
    {
      psrc = cols [c].start ;

      /* move and compact the column */
      _assert (pdest <= psrc) ;
      cols [c].start = pdest - 0 ;
      length = cols [c].length ;
      for (j = 0 ; j < length ; j++)
      {
        r = A [psrc++] ;
        if (_rowIsAlive (rows, r))
        {
          A [pdest] = r ;
        }
      }
      cols [c].length = pdest - cols [c].start ;
    }
  }

  /* === Prepare to defragment the rows =============================== */

  for (r = 0 ; r < n_row ; r++)
  {
    if (_rowIsDead (rows, r) || (rows [r].length == 0))
    {
      /* This row is already dead, or is of zero length.  Cannot compact
       * a row of zero length, so kill it.  NOTE: in the current version,
       * there are no zero-length live rows.  Kill the row (for the first
       * time, or again) just to be safe. */
      _killRow (rows, r) ;
    }
    else
    {
      /* save first column index in Row [r].shared2.first_column */
      psrc = rows [r].start ;
      rows [r].first_column = A [psrc] ;
      _assert (_rowIsAlive (rows, r)) ;
      /* flag the start of the row with the one's complement of row */
      A [psrc] = _onesComplement (r) ;
      if (!_ndebug)
      {
        debug_rows++ ;
      } /* NDEBUG */
    }
  }

  /* === Defragment the rows ========================================== */

  psrc = pdest ;
  while (psrc < pfree)
  {
    /* find a negative number ... the start of a row */
    if (A [psrc++] < 0)
    {
      psrc-- ;
      /* get the row index */
      r = _onesComplement (A [psrc]) ;
      _assert (r >= 0 && r < n_row) ;
      /* restore first column index */
      A [psrc] = rows [r].first_column ;
      _assert (_rowIsAlive (rows, r)) ;
      _assert (rows [r].length > 0) ;
      /* move and compact the row */
      _assert (pdest <= psrc) ;
      rows [r].start = pdest - 0 ;
      length = rows [r].length ;
      for (j = 0 ; j < length ; j++)
      {
        c = A [psrc++] ;
        if (_colIsAlive (cols, c))
        {
          A [pdest++] = c ;
        }
      }
      rows [r].length = pdest - rows [r].start ;
      _assert (rows [r].length > 0) ;
      if (!_ndebug)
      {
        debug_rows-- ;
      } /* NDEBUG */
    }
  }
  /* ensure we found all the rows */
  _assert (debug_rows == 0) ;

  /* === Return the new value of pfree ================================ */

  return (pdest - 0) ;
}


/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

/// Clears the Row [].shared2.mark array, and returns the new tag_mark.
///
/// @param tag_mark new value of tag_mark
/// @param max_mark max allowed value of tag_mark
/// @param n_row number of rows in A
/// @param Row Row [0 ... n_row-1].shared2.mark is set to zero
/// @return the new value for tag_mark
int _clearMark (int tag_mark, int max_mark, int n_row,
    List<Row> Row)
{
  if (tag_mark <= 0 || tag_mark >= max_mark)
  {
    for (int r = 0 ; r < n_row ; r++)
    {
      if (_rowIsAlive (Row, r))
      {
        Row [r].mark = 0 ;
      }
    }
    tag_mark = 1 ;
  }

  return (tag_mark) ;
}


/* ========================================================================== */
/* === print_report ========================================================= */
/* ========================================================================== */

void _printReport(String method, List<int> stats)
{
  _print ("\n$method version $MAIN_VERSION.$SUB_VERSION, $DATE: ") ;

  if (stats == null)
  {
    _print ("No statistics available.\n") ;
    return ;
  }

  final i1 = stats [INFO1] ;
  final i2 = stats [INFO2] ;
  final i3 = stats [INFO3] ;

  if (stats [STATUS] >= 0)
  {
    _print ("OK.  ") ;
  }
  else
  {
    _print ("ERROR.  ") ;
  }

  switch (stats [STATUS])
  {

    case OK_BUT_JUMBLED:

      _print("Matrix has unsorted or duplicate row indices.\n") ;

      _print("$method: number of duplicate or out-of-order row indices: $i3\n") ;

      _print("$method: last seen duplicate or out-of-order row index:   ${_index (i2)}\n") ;

      _print("$method: last seen in column:                             ${_index (i1)}") ;


      _print("\n") ;

      _print("$method: number of dense or empty rows ignored:           ${stats [DENSE_ROW]}\n") ;

      _print("$method: number of dense or empty columns ignored:        ${stats [DENSE_COL]}\n") ;

      _print("$method: number of garbage collections performed:         ${stats [DEFRAG_COUNT]}\n") ;
      break ;

    case OK:

      _print("\n") ;

      _print("$method: number of dense or empty rows ignored:           ${stats [DENSE_ROW]}\n") ;

      _print("$method: number of dense or empty columns ignored:        ${stats [DENSE_COL]}\n") ;

      _print("$method: number of garbage collections performed:         ${stats [DEFRAG_COUNT]}\n") ;
      break ;

    case ERROR_A_NOT_PRESENT:

      _print("Array A (row indices of matrix) not present.\n") ;
      break ;

    case ERROR_P_NOT_PRESENT:

      _print("Array p (column pointers for matrix) not present.\n") ;
      break ;

    case ERROR_NROW_NEGATIVE:

      _print("Invalid number of rows ($i1).\n") ;
      break ;

    case ERROR_NCOL_NEGATIVE:

      _print("Invalid number of columns ($i1).\n") ;
      break ;

    case ERROR_NNZ_NEGATIVE:

      _print("Invalid number of nonzero entries ($i1).\n") ;
      break ;

    case ERROR_P0_NONZERO:

      _print("Invalid column pointer, p [0] = $i1, must be zero.\n");
      break ;

    case ERROR_A_TOO_SMALL:

      _print("Array A too small.\n") ;
      _print("        Need Alen >= $i1, but given only Alen = $i2.\n") ;
      break ;

    case ERROR_COL_LENGTH_NEGATIVE:

      _print
      ("Column ${_index (i1)} has a negative number of nonzero entries ($i2).\n") ;
      break ;

    case ERROR_ROW_INDEX_OUT_OF_BOUNDS:

      _print
      ("Row index (row ${_index (i2)}) out of bounds (${_index (0)} to ${_index (i3-1)}) in column ${_index (i1)}.\n") ;
      break ;

    case ERROR_OUT_OF_MEMORY:

      _print("Out of memory.\n") ;
      break ;

    /* v2.4: internal-error case deleted */
  }
}


/* ========================================================================== */
/* === colamd debugging routines ============================================ */
/* ========================================================================== */

/// At this point, all empty rows and columns are dead.  All live columns
/// are "clean" (containing no dead rows) and simplicial (no supercolumns
/// yet).  Rows may contain dead columns, but all live rows contain at
/// least one live column.
void _debugStructures (int n_row, int n_col,
    List<Row> rows, List<Col> cols, Int32List A, int n_col2)
{
  if (!_ndebug)
  {
    /* === Local variables ============================================== */

    int i, c, cp, cp_end, len, score, r, rp, rp_end, deg ;

    /* === Check A, Row, and Col ======================================== */

    for (c = 0 ; c < n_col ; c++)
    {
      if (_colIsAlive (cols, c))
      {
        len = cols [c].length ;
        score = cols [c].score ;
        _debug4 ("initial live col $c $len $score\n") ;
        _assert (len > 0) ;
        _assert (score >= 0) ;
        _assert (cols [c].thickness == 1) ;
        cp = cols [c].start ;
        cp_end = cp + len ;
        while (cp < cp_end)
        {
          r = A [cp++] ;
          _assert (_rowIsAlive (rows, r)) ;
        }
      }
      else
      {
        i = cols [c].order ;
        _assert (i >= n_col2 && i < n_col) ;
      }
    }

    for (r = 0 ; r < n_row ; r++)
    {
      if (_rowIsAlive (rows, r))
      {
        i = 0 ;
        len = rows [r].length ;
        deg = rows [r].degree ;
        _assert (len > 0) ;
        _assert (deg > 0) ;
        rp = rows [r].start ;
        rp_end = rp + len ;
        while (rp < rp_end)
        {
          c = A [rp++] ;
          if (_colIsAlive (cols, c))
          {
            i++ ;
          }
        }
        _assert (i > 0) ;
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
void _debugDegLists (int n_row, int n_col,
    List<Row> rows, List<Col> cols, Int32List head, int min_score,
    int should, int max_deg)
{
  if (!_ndebug)
  {
    /* === Local variables ============================================== */

    int deg, col, have, row ;

    /* === Check the degree lists ======================================= */

    if (n_col > 10000 && debugLevel <= 0)
    {
      return ;
    }
    have = 0 ;
    _debug4 ("Degree lists: $min_score\n") ;
    for (deg = 0 ; deg <= n_col ; deg++)
    {
      col = head [deg] ;
      if (col == _empty)
      {
        continue ;
      }
      _debug4 ("$deg:") ;
      while (col != _empty)
      {
        _debug4 (" $col") ;
        have += cols [col].thickness ;
        _assert (_colIsAlive (cols, col)) ;
        col = cols [col].degree_next ;
      }
      _debug4 ("\n") ;
    }
    _debug4 ("should $should have $have\n") ;
    _assert (should == have) ;

    /* === Check the row degrees ======================================== */

    if (n_row > 10000 && debugLevel <= 0)
    {
      return ;
    }
    for (row = 0 ; row < n_row ; row++)
    {
      if (_rowIsAlive (rows, row))
      {
        _assert (rows [row].degree <= max_deg) ;
      }
    }
  }
}


/* ========================================================================== */
/* === debug_mark =========================================================== */
/* ========================================================================== */

/// Ensures that the tag_mark is less that the maximum and also ensures that
/// each entry in the mark array is less than the tag mark.
void _debugMark (int n_row, List<Row> rows, int tag_mark,
    int max_mark)
{
  if (!_ndebug)
  {
    /* === Check the Row marks ========================================== */

    _assert (tag_mark > 0 && tag_mark <= max_mark) ;
    if (n_row > 10000 && debugLevel <= 0)
    {
      return ;
    }
    for (int r = 0 ; r < n_row ; r++)
    {
      _assert (rows [r].mark < tag_mark) ;
    }
  }
}


/* ========================================================================== */
/* === debug_matrix ========================================================= */
/* ========================================================================== */

/// Prints out the contents of the columns and the rows.
void _debugMatrix (int n_row, int n_col,
    List<Row> rows, List<Col> cols, List<int> A)
{
  if (!_ndebug)
  {
    /* === Local variables ============================================== */

    int r ;
    int c ;
    int rp ;
    int rp_end ;
    int cp ;
    int cp_end ;

    /* === Dump the rows and columns of the matrix ====================== */

    if (debugLevel < 3)
    {
      return ;
    }
    _debug3 ("DUMP MATRIX:\n") ;
    for (r = 0 ; r < n_row ; r++)
    {
      _debug3 ("Row $r alive? ${_rowIsAlive (rows, r) ? 1 : 0}\n") ;
      if (_rowIsDead (rows, r))
      {
        continue ;
      }
      _debug3 ("start ${rows [r].start} length ${rows [r].length} degree ${rows [r].degree}\n") ;
      rp = rows [r].start ;
      rp_end = rp + rows [r].length ;
      while (rp < rp_end)
      {
        c = A [rp++] ;
        _debug4 ("  ${_colIsAlive (cols, c) ? 1 : 0} col $c\n") ;
      }
    }

    for (c = 0 ; c < n_col ; c++)
    {
      _debug3 ("Col $c alive? ${_colIsAlive (cols, c) ? 1 : 0}\n") ;
      if (_colIsDead (cols, c))
      {
        continue ;
      }
      _debug3 ("start ${cols [c].start} length ${cols [c].length} shared1 ${cols [c].thickness} shared2 ${cols [c].score}\n") ;
      cp = cols [c].start ;
      cp_end = cp + cols [c].length ;
      while (cp < cp_end)
      {
        r = A [cp++] ;
        _debug4 ("  ${_rowIsAlive (rows, r) ? 1 : 0} row $r\n") ;
      }
    }
  }
}
