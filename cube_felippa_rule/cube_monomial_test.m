function cube_monomial_test ( degree_max )

%*****************************************************************************80
%
%% CUBE_MONOMIAL_TEST tests CUBE_MONOMIAL.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    04 September 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer DEGREE_MAX, the maximum total degree of the
%    monomials to check.
%
  a(1:3) = -1.0;
  b(1:3) = +1.0;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'CUBE_MONOMIAL_TEST\n' );
  fprintf ( 1, '  For a cube,\n' );
  fprintf ( 1, '  CUBE_MONOMIAL returns the exact value of the\n' );
  fprintf ( 1, '  integral of X^ALPHA Y^BETA Z^GAMMA\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Volume = %f\n', cube_volume ( a, b ) );
  fprintf ( 1, '\n' );
  fprintf ( 1, '     ALPHA      BETA     GAMMA      INTEGRAL\n' );
  fprintf ( 1, '\n' );

  for alpha = 0 : degree_max
    expon(1) = alpha;
    for beta = 0 : degree_max - alpha
      expon(2) = beta;
      for gamma = 0 : degree_max - alpha - beta
        expon(3) = gamma;
        value = cube_monomial ( a, b, expon );
        fprintf ( 1, '  %8d  %8d  %8d  %14.6e\n', expon(1:3), value );
      end
    end
  end

  return
end
