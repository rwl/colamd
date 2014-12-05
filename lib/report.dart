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
int debugLevel = 0;

void _debug0(String format) {
  if (!_ndebug) {
    _print(format);
  }
}

void _debug1(String format) {
  if (!_ndebug) {
    if (debugLevel >= 1) _print(format);
  }
}

void _debug2(String format) {
  if (!_ndebug) {
    if (debugLevel >= 2) _print(format);
  }
}

void _debug3(String format) {
  if (!_ndebug) {
    if (debugLevel >= 3) _print(format);
  }
}

void _debug4(String format) {
  if (!_ndebug) {
    if (debugLevel >= 4) _print(format);
  }
}

void _assert(bool a) {
  if (!_ndebug) {
    assert(a);
  }
}

void _assertInt(int a) {
  _assert(a != 0);
}

void colamdReport(List<int> stats) {
  _printReport("colamd", stats);
}

void symamdReport(List<int> stats) {
  _printReport("symamd", stats);
}

void _printReport(String method, List<int> stats) {
  var s = "\n$method version $MAIN_VERSION.$SUB_VERSION, $DATE: ";

  if (stats == null) {
    _print(s + "No statistics available.");
    return;
  }

  final i1 = stats[INFO1];
  final i2 = stats[INFO2];
  final i3 = stats[INFO3];

  if (stats[STATUS] >= 0) {
    _print("$s OK.  ");
  } else {
    _print("$s ERROR.  ");
  }

  switch (stats[STATUS]) {

    case OK_BUT_JUMBLED:

      _print("Matrix has unsorted or duplicate row indices.");

      _print("$method: number of duplicate or out-of-order row indices: $i3");

      _print("$method: last seen duplicate or out-of-order row index:   ${_index (i2)}");

      _print("$method: last seen in column:                             ${_index (i1)}");


      _print("");

      _print("$method: number of dense or empty rows ignored:           ${stats [DENSE_ROW]}");

      _print("$method: number of dense or empty columns ignored:        ${stats [DENSE_COL]}");

      _print("$method: number of garbage collections performed:         ${stats [DEFRAG_COUNT]}");
      break;

    case OK:

      _print("");

      _print("$method: number of dense or empty rows ignored:           ${stats [DENSE_ROW]}");

      _print("$method: number of dense or empty columns ignored:        ${stats [DENSE_COL]}");

      _print("$method: number of garbage collections performed:         ${stats [DEFRAG_COUNT]}");
      break;

    case ERROR_A_NOT_PRESENT:

      _print("Array A (row indices of matrix) not present.");
      break;

    case ERROR_P_NOT_PRESENT:

      _print("Array p (column pointers for matrix) not present.");
      break;

    case ERROR_NROW_NEGATIVE:

      _print("Invalid number of rows ($i1).");
      break;

    case ERROR_NCOL_NEGATIVE:

      _print("Invalid number of columns ($i1).");
      break;

    case ERROR_NNZ_NEGATIVE:

      _print("Invalid number of nonzero entries ($i1).");
      break;

    case ERROR_P0_NONZERO:

      _print("Invalid column pointer, p [0] = $i1, must be zero.");
      break;

    case ERROR_A_TOO_SMALL:

      _print("Array A too small.");
      _print("        Need Alen >= $i1, but given only Alen = $i2.");
      break;

    case ERROR_COL_LENGTH_NEGATIVE:

      _print("Column ${_index (i1)} has a negative number of nonzero entries ($i2).");
      break;

    case ERROR_ROW_INDEX_OUT_OF_BOUNDS:

      _print("Row index (row ${_index (i2)}) out of bounds (${_index (0)} to ${_index (i3-1)}) in column ${_index (i1)}.");
      break;

    case ERROR_OUT_OF_MEMORY:

      _print("Out of memory.");
      break;

  }
}
