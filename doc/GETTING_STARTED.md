# Getting Started

If you have reached this documentation page, it is assumed that you are interested in digging into
the details of the implementation. Otherwise, please refer to the
<a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/README.md">Root README</a> file
for a basic overview.

### Introduction to the Implementation

The best way to introduce yourself to the code is by consulting the unit and integration tests found
in sub-directories of the
<a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/src/testing">testing</a>
directory. Unit test can be used to understand isolated functionality while the available
integration tests build from verification of basic preprocessor functionality to full-scale
convergence order testing of the implemented methods. The full list of test functions can be found
in the <a href="test.html">Test List</a>.

Consulting the following tests in the order specified below should result in a gradual exposure to
the code features:
- <a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/src/testing/integration/test_integration_mesh.c">test_integration_mesh</a>;
- <a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/src/testing/integration/test_integration_fe_init.c">test_integration_fe_init</a>;
- <a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/src/testing/integration/test_integration_geometry.c">test_integration_geometry</a>;
- <a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/src/testing/integration/test_integration_linearization.c">test_integration_linearization</a>;
- <a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/src/testing/integration/test_integration_convergence.c">test_integration_convergence</a>.

Note that many tests have associated `.data` files in the
<a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master/input/testing">testing input directory</a>
which are compared with data structures generated in the code. Of particular note is that these data
files are easily accessible to human readers, further facilitating understanding of code
functionality.


### Adding Your Own Contributions

It is highly encouraged to add any of your own contributions through a *test-driven development*
process by first implementing a test for your desired feature and then writing the associated code
needed to pass the test. In this manner, your knowledge and expertise will be embedded directly in
the code through the test itself and the associated documentation as you perform the implementation.
