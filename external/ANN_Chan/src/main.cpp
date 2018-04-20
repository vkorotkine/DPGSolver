// Timothy Chan 12/05
// approximate nearest neighbors: the SSS method (static version)

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <limits>
#define sq(x) (((float) (x))*((float) (x)))
#define MAX (1<<29)
#define FLT_MAX std::numeric_limits<float>::max()

using namespace std;

typedef int* Point;
int d, shift;
float eps, r, r_sq;
Point ans, q1, q2;

inline int less_msb(int x, int y) { return x < y && x < (x^y); }

int cmp_shuffle(Point* p, Point* q) {
	int j, k, x, y;
	for (j = k = x = 0; k < d; k++)
		if (less_msb(x, y = ((*p)[k]+shift)^((*q)[k]+shift))) {
			j = k; x = y;
		}
	return (*p)[j]-(*q)[j];
}

void SSS_preprocess(Point P[], int n) {
	shift = (int) (drand48()*MAX);
	q1 = new int[d]; q2 = new int[d];
	qsort((void*) P, n, sizeof(Point),(int (*)(const void *, const void *)) cmp_shuffle);
}

void check_dist(Point p, Point q)
{
	int j; float z;
	for (j = 0, z = 0; j < d; j++) z += sq(p[j]-q[j]);
	if (z < r_sq) {
		r_sq = z; r = sqrt(z); ans = p;
		for (j = 0; j < d; j++) {
			q1[j] = (q[j]>r) ? (q[j]-(int)ceil(r)) : 0;
			q2[j] = (q[j]+r<MAX) ? (q[j]+(int)ceil(r)) : MAX;
		}
	}
}

float dist_sq_to_box(Point q, Point p1, Point p2) {
	int i, j, x, y; float z;
	for (j = x  = 0; j < d; j++)
		if (less_msb(x,y = (p1[j]+shift)^(p2[j]+shift))) x = y;
	frexp(x,&i);
	for (j = 0, z = 0; j < d; j++) {
		x = ((p1[j]+shift)>>1)<<i; y = x+(1<<i);
		if (q[j]+shift < x) z += sq(q[j]+shift-x);
		else if (q[j]+shift > y) z += sq(q[j]+shift-y);
	}
	return z;
}

void SSS_query0(Point P[], int n, Point q) {
	if (n == 0) return;
	check_dist(P[n/2],q);
	if (n == 1 || dist_sq_to_box(q, P[0],P[n-1])*sq(1+eps) > r_sq) return;
	if (cmp_shuffle(&q, &P[n/2]) < 0) {
		SSS_query0(P,n/2,q);
		if (cmp_shuffle(&q2, &P[n/2]) > 0) SSS_query0(P+n/2+1, n-n/2-1, q);
	} else {
		SSS_query0(P+n/2+1, n-n/2-1, q);
		if (cmp_shuffle(&q1,&P[n/2]) < 0) SSS_query0(P,n/2,q);
	}
}

Point SSS_query(Point* P, int n, Point q) {
	r_sq = FLT_MAX;
	SSS_query0(P,n,q);
	return ans;
}

#define n 3
#define m 2
#define d 2

int main (int argc, char** argv)
{
	const float eps = 1e-4;

	const int p_data[n][d] = 
		{ { 1, 2, }, {4, 3,}, {10, 5,},};
	const int q_data[m][d] = 
		{ { 0, 0, }, {4, 4,}, };
	
	srand48(12121+n+m+d);
	Point* P;
	P = new Point[n];
	for (int i = 0, ind = 0; i < n; i++) {
		P[i] = new int[d];
		for (int j = 0; j < d; j++)
			P[i][j] = p_data[i][ind++];
	}
	SSS_preprocess(P,n);

	Point q;
	q = new int[d];
	for (int i = 0, ind = 0; i < m; i++) {
		for (int j = 0; j < d; j++)
			q[j] = q_data[i][ind++];
		SSS_query(P,n,q);
		cout << r << r_sq << "\n";
	}
	for (int i = 0; i < n; i++) delete P[i];
	delete P; delete q;
}

/*

Algorithm to sort the points.

ImmutableList<Point> OrderByDistance(Point start, ImmutableSet<Point> points)
{
  var current = start;
  var remaining = points;
  var path = ImmutableList<Point>.Empty.Add(start);
  while(!remaining.IsEmpty)
  {
    var next = Closest(current, remaining);
    path = path.Add(next);
    remaining = remaining.Remove(next);
    current = next;
  }
  return path;
}
*/
