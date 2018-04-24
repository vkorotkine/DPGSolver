// Timothy Chan 12/05
// approximate nearest neighbors: the SSS method (static version)

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <limits>
#define sq(x) (((float) (x))*((float) (x)))
//#define MAX (1<<29)
#define MAX INT_MAX
#define FLT_MAX std::numeric_limits<float>::max()
#define EXIT ({ fflush(stdout); abort(); })

using namespace std;

#define d 2

typedef int* Point;
int shift;
float eps, r, r_sq;
Point ans, q1, q2;

// p.4: Used to avoid computing msb: msb(x) < msb(y) iff ( x < y && x < (x^y) ). Cites ref [5].
inline int less_msb(int x, int y) { return x < y && x < (x^y); }

/** \brief Comparison function for Points.
 *
 *  p <~ q iff x_j <= y_j where j is the index of the 'm'ost 's'ignificant 'b'it of (x xor y) for j \in d. The <~
 *  operator, called the shuffle order, is used to denote comparison of the shuffle of two points written in binary.
 */
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
	// l2 norm diff
	for (j = 0, z = 0; j < d; j++) {
//printf("%d %d %d\n",p[j],q[j],p[j]-q[j]);
		z += sq(p[j]-q[j]);
	}
//printf("% .3e\n",z);

	if (z < r_sq) {
		r_sq = z; r = sqrt(z); ans = p; // p.3 line 2
		for (j = 0; j < d; j++) {
			q1[j] = (q[j]>r) ? (q[j]-(int)ceil(r)) : 0;       // q1 == q^{s-[r]} (p.3, line 9)
			q2[j] = (q[j]+r<MAX) ? (q[j]+(int)ceil(r)) : MAX; // q2 == q^{s+[r]} (p.3, line 6)
		}
	}
}

float dist_sq_to_box(Point q, Point p1, Point p2) {
	int i, j, x, y; float z;
	for (j = x  = 0; j < d; j++)
		if (less_msb(x,y = (p1[j]+shift)^(p2[j]+shift))) x = y;
	frexp(x,&i); // Extract most significant bit
	for (j = 0, z = 0; j < d; j++) {
		x = ((p1[j]+shift)>>i)<<i; y = x+(1<<i);
		if (q[j]+shift < x) z += sq(q[j]+shift-x);
		else if (q[j]+shift > y) z += sq(q[j]+shift-y);
	}
	return z;
}

void SSS_query0(Point P[], int n, Point q)
{
//static int count = 0;
//printf("n: %d %d %d %d\n",count++,n,P[0][0],P[0][1]);
	if (n == 0)
		return;
	check_dist(P[n/2],q); // p.3 line 1 (a = 0, b = n => (a+b)/2 = n/2)

	printf("% 3d %d %d %d %d % .3e % .3e\n",
	       n,P[0][0],P[0][1],P[n-1][0],P[n-1][1],dist_sq_to_box(q, P[0],P[n-1])*sq(1+eps),r_sq);

	if (n == 1 || dist_sq_to_box(q, P[0],P[n-1])*sq(1+eps) > r_sq) return; // p.3 line 3
	if (cmp_shuffle(&q, &P[n/2]) < 0) { // p.3 line 4
//printf("lower\n");
		SSS_query0(P,n/2,q);          // p.3 line 5 (binary search in lower half)
		if (cmp_shuffle(&q2, &P[n/2]) > 0) SSS_query0(P+n/2+1, n-n/2-1, q); // p.3 line 6
	} else {
//printf("higher\n");
		SSS_query0(P+n/2+1, n-n/2-1, q); // p.3 line 8 (binary search in upper half)
//printf("h2: %d %d %d %d %d\n",q1[0],q1[1],P[n/2][0],P[n/2][1],
//       cmp_shuffle(&q1,&P[n/2]) < 0);
		if (cmp_shuffle(&q1,&P[n/2]) < 0) SSS_query0(P,n/2,q); // p.3 line 9
	}
}

Point SSS_query(Point* P, int n, Point q) {
	r_sq = FLT_MAX;
	SSS_query0(P,n,q);
	return ans;
}

#define n 16
#define m 3

int main (int argc, char** argv)
{
	const float eps = 0.0;

	const int p_scale = 4;
	const int p_data[n][d] =
		{ {0,0,}, {0,1,}, {0,2,}, {0,3,},
		  {1,0,}, {1,1,}, {1,2,}, {1,3,},
		  {2,0,}, {2,1,}, {2,2,}, {2,3,},
		  {3,0,}, {3,1,}, {3,2,}, {3,3,},
	      };
	const int q_data[m][d] =
		{ {13,9,}, {4,5,}, {0,8,}, };

	srand48(31415+n+m+d);
	Point* P;
	P = new Point[n];
	for (int i = 0; i < n; i++) {
		P[i] = new int[d];
		for (int j = 0; j < d; j++)
			P[i][j] = p_scale*p_data[i][j];
	}

	SSS_preprocess(P,n);

	Point q;
	q = new int[d];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < d; j++)
			q[j] = q_data[i][j];
		SSS_query(P,n,q);
		cout << "r, r2, ans: " << r << ", " << r_sq << ", (" << ans[0] << ","<< ans[1] << ")" << "\n";
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
