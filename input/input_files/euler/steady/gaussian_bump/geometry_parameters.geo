/** Geometry parameters for test case: euler/steady/gaussian_bump
 *
 *  Note: As x_max was increased, problems with convergence were encountered 
 *        when too few elements were used to represent the curved boundary 
 *        geometry (Seemingly related to the discrete curvature constraint).
 */

a = 0.2;     // height (0.4 was giving good results)
b = -0.0314; // x-offset
c = 0.5;     // width
d = 0.1;     // x-linear scaling

h = 2.0;
x_max = 3.0;
