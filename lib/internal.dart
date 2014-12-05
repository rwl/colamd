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

num _max(num a, num b) => a > b ? a : b;

num _min(num a, num b) => a < b ? a : b;

int _denseDegree(double alpha, int n) {
  return _max(16.0, (alpha) * math.sqrt(n.toDouble())).toInt();
}

int _onesComplement(int r) => -(r) - 1;

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
bool _rowIsDead(List<Row> Row, int r) {
  return _rowIsMarkedDead(Row[r].mark);
}

bool _rowIsMarkedDead(int row_mark) {
  return (row_mark < _alive);
}

bool _rowIsAlive(List<Row> Row, int r) {
  return (Row[r].mark >= _alive);
}

bool _colIsDead(List<Col> Col, int c) {
  return (Col[c].start < _alive);
}

bool _colIsAlive(List<Col> Col, int c) {
  return (Col[c].start >= _alive);
}

bool _colIsDeadPrincipal(List<Col> Col, int c) {
  return (Col[c].start == deadPrincipal);
}

void _killRow(List<Row> Row, int r) {
  Row[r].mark = _dead;
}

void _killPrincipalCol(List<Col> Col, int c) {
  Col[c].start = deadPrincipal;
}

void _killNonPrincipalCol(List<Col> Col, int c) {
  Col[c].start = deadNonPrincipal;
}
