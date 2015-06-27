function bohach2_test ( )

%*****************************************************************************80
%
%% BOHACH2_TEST works with the Bohachevsky function #2.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 December 2011
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Zbigniew Michalewicz,
%    Genetic Algorithms + Data Structures = Evolution Programs,
%    Third Edition,
%    Springer Verlag, 1996,
%    ISBN: 3-540-60676-9,
%    LC: QA76.618.M53.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'BOHACH2_TEST:\n' );
  fprintf ( 1, '  Test COORDINATE_SEARCH with the Bohachevsky function #2.\n' );
  n = 2;

  x = [ 0.6, 1.3 ];
  r8vec_print ( n, x, '  Initial point X0:' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  F(X0) = %g\n', bohach2 ( x ) );

  flag = 0;
  x = coordinate_search ( x, @bohach2, flag );
  r8vec_print ( n, x, '  Estimated minimizer X1:' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  F(X1) = %g\n', bohach2 ( x ) );
%
%  Demonstrate correct minimizer.
%
  x = [ 0.0, 0.0 ];
  r8vec_print ( n, x, '  Correct minimizer X*:' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  F(X*) = %g\n', bohach2 ( x ) );

  return
end
