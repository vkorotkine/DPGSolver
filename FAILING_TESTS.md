# Documented Failing Tests

Tests which are failing either for an expected reason or because further investigation is required are listed here. Note
that the tests were run in the release build. All tests not documented here should be passing; please report an issue if
this is not the case.
### test_integration_restart___integration/restart/TEST_Euler_SupersonicVortex_DG_ParametricTri2D

Please see the comments in [test_integration_restart.c](src/testing/integration/test_integration_restart.c).

### test_integration_convergence___euler/joukowski/TEST_Euler_Joukowski_DG_ParametricMixed2D

The convergence orders are tending towards the expected optimal orders but require finer meshes than those
currently used to be within the required tolerance as specified by the 'ACCEPTABLE_DISCOUNT' parameter in
[test_integration_convergence_support.c](src/testing/integration/test_integration_convergence_support.c).
