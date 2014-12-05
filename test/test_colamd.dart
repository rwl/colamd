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
library colamd.test;

import 'dart:io' show stdout;
import 'dart:typed_data';
import 'package:unittest/unittest.dart';
import 'package:colamd/colamd.dart' as colamd;

const int A_NNZ = 11;
const int A_NROW = 5;
const int A_NCOL = 4;
const int ALEN = 150;

const int B_NNZ = 4;
const int B_N = 5;

/// COLAMD / SYMAMD example
///
/// colamd example of use, to order the columns of a 5-by-4 matrix with
/// 11 nonzero entries in the following nonzero pattern, with default knobs.
///
///   x 0 x 0
///   x 0 x x
///   0 x x 0
///   0 0 x x
///   x x 0 0
///
/// symamd example of use, to order the rows and columns of a 5-by-5
/// matrix with 13 nonzero entries in the following nonzero pattern,
/// with default knobs.
///
///   x x 0 0 0
///   x x x x 0
///   0 x x 0 0
///   0 x 0 x x
///   0 0 0 x x
///
/// (where x denotes a nonzero value).
main() {

  colamd.debug = true;
  colamd.nprint = false;
  colamd.debugLevel = 0;

  test('colamd', () {

    /* ====================================================================== */
    /* input matrix A definition */
    /* ====================================================================== */

    final A = new Int32List(ALEN);
    final AA = [0, 1, 4, // row indices of nonzeros in column 0
      2, 4, // row indices of nonzeros in column 1
      0, 1, 2, 3, // row indices of nonzeros in column 2
      1, 3]; // row indices of nonzeros in column 3
    //System.arraycopy(AA, 0, A, 0, AA.length) ;
    A.setAll(0, AA);

    final p = [0, // column 0 is in A [0..2]
      3, // column 1 is in A [3..4]
      5, // column 2 is in A [5..8]
      9, // column 3 is in A [9..10]
      A_NNZ]; // number of nonzeros in A

    /* ====================================================================== */
    /* input matrix B definition */
    /* ====================================================================== */

    final B = [ // Note: only strictly lower triangular part
      // is included, since symamd ignores the
      // diagonal and upper triangular part of B.

      1, // row indices of nonzeros in column 0
      2, 3, // row indices of nonzeros in column 1
      // row indices of nonzeros in column 2 (none)
      4 // row indices of nonzeros in column 3
    ]; // row indices of nonzeros in column 4 (none)

    final q = [0, // column 0 is in B [0]
      1, // column 1 is in B [1..2]
      3, // column 2 is empty
      3, // column 3 is in B [3]
      4, // column 4 is empty
      B_NNZ]; // number of nonzeros in strictly lower B

    /* ====================================================================== */
    /* other variable definitions */
    /* ====================================================================== */

    final perm = new Int32List(B_N + 1); // note the size is N+1
    final stats = new Int32List(colamd.STATS); // for colamd and symamd output statistics

    int row, col, pp, length, ok;

    /* ====================================================================== */
    /* dump the input matrix A */
    /* ====================================================================== */

    stdout.write("colamd $A_NROW-by-$A_NCOL input matrix:\n");
    for (col = 0; col < A_NCOL; col++) {
      length = p[col + 1] - p[col];
      stdout.write("Column $col, with $length entries:\n");
      for (pp = p[col]; pp < p[col + 1]; pp++) {
        row = A[pp];
        stdout.write("    row $row\n");
      }
    }

    /* ====================================================================== */
    /* order the matrix.  Note that this destroys A and overwrites p */
    /* ====================================================================== */

    ok = colamd.colamd(A_NROW, A_NCOL, ALEN, A, p, null, stats);
    colamd.colamdReport(stats);

    if (ok == 0) {
      stdout.write("colamd error!\n");
      fail('');
    }

    /* ====================================================================== */
    /* stdout.write the column ordering */
    /* ====================================================================== */

    stdout.write("colamd column ordering:\n");
    stdout.write("1st column: ${p [0]}\n");
    stdout.write("2nd column: ${p [1]}\n");
    stdout.write("3rd column: ${p [2]}\n");
    stdout.write("4th column: ${p [3]}\n");

    expect(1, equals(p[0]));
    expect(0, equals(p[1]));
    expect(2, equals(p[2]));
    expect(3, equals(p[3]));

    /* ====================================================================== */
    /* dump the strictly lower triangular part of symmetric input matrix B */
    /* ====================================================================== */

    stdout.write("\n\nsymamd $B_N-by-$B_N input matrix:\n");
    stdout.write("Entries in strictly lower triangular part:\n");
    for (col = 0; col < B_N; col++) {
      length = q[col + 1] - q[col];
      stdout.write("Column $col, with $length entries:\n");
      for (pp = q[col]; pp < q[col + 1]; pp++) {
        row = B[pp];
        stdout.write("    row $row\n");
      }
    }

    /* ====================================================================== */
    /* order the matrix B.  Note that this does not modify B or q. */
    /* ====================================================================== */

    ok = colamd.symamd(B_N, B, q, perm, null, stats);
    colamd.symamdReport(stats);

    if (ok == 0) {
      stdout.write("symamd error!\n");
      fail('');
    }

    /* ====================================================================== */
    /* stdout.write the symmetric ordering */
    /* ====================================================================== */

    stdout.write("symamd column ordering:\n");
    stdout.write("1st row/column: ${perm [0]}\n");
    stdout.write("2nd row/column: ${perm [1]}\n");
    stdout.write("3rd row/column: ${perm [2]}\n");
    stdout.write("4th row/column: ${perm [3]}\n");
    stdout.write("5th row/column: ${perm [4]}\n");

    expect(0, equals(perm[0]));
    expect(2, equals(perm[1]));
    expect(1, equals(perm[2]));
    expect(3, equals(perm[3]));
    expect(4, equals(perm[4]));
  });
}
