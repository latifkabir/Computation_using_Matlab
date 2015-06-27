function fem1d_bvp_quadratic_test08 ( )

%*****************************************************************************80
%
%% FEM1D_BVP_QUADRATIC_TEST08 carries out test case #8.
%
%  Discussion:
%
%    Use A8, C8, F8, EXACT8, EXACT_UX8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 June 2014
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Dianne O'Leary,
%    Scientific Computing with Case Studies,
%    SIAM, 2008,
%    ISBN13: 978-0-898716-66-5,
%    LC: QA401.O44.
%
  n = 11;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'FEM1D_BVP_QUADRATIC_TEST08\n' );
  fprintf ( 1, '  Solve -( A(x) U''(x) )'' + C(x) U(x) = F(x)\n' );
  fprintf ( 1, '  for 0 < x < 1, with U(0) = U(1) = 0.\n' );
  fprintf ( 1, '  A8(X)  = 1.0\n' );
  fprintf ( 1, '  C8(X)  = 0.0\n' );
  fprintf ( 1, '  F8(X)  = X * ( X + 3 ) * exp ( X ),   X <= 2/3\n' );
  fprintf ( 1, '         = 2 * exp ( 2/3),                   2/3 < X\n' );
  fprintf ( 1, '  U8(X)  = X * ( 1 - X ) * exp ( X ),   X <= 2/3\n' );
  fprintf ( 1, '         = X * ( 1 - X ) * exp ( 2/3 ),      2/3 < X\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of nodes = %d\n', n );
%
%  Geometry definitions.
%
  x_first = 0.0;
  x_last = 1.0;
  x = linspace ( x_first, x_last, n );
  x = x(:);

  u = fem1d_bvp_quadratic ( n, @a8, @c8, @f8, x );

  uexact = exact8 ( x );

  fprintf ( 1, '\n' );
  fprintf ( 1, '     I    X         U         Uexact    Error\n' );
  fprintf ( 1, '\n' );

  for i = 1 : n
    fprintf ( 1, '  %4d  %8f  %8f  %8f  %8e\n', ...
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) ) );
  end

  e1 = l1_error ( n, x, u, @exact8 );
  e2 = l2_error_quadratic ( n, x, u, @exact8 );
  h1s = h1s_error_quadratic ( n, x, u, @exact_ux8 );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  l1 norm of error  = %g\n', e1 );
  fprintf ( 1, '  L2 norm of error  = %g\n', e2 );
  fprintf ( 1, '  Seminorm of error = %g\n', h1s );

  return
end
function value = a8 ( x )

%*****************************************************************************80
%
%% A8 evaluates A function #8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 June 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the evaluation point.
%
%    Output, real VALUE, the value of A(X).
%
  value = 1.0;

  return
end
function value = c8 ( x )

%*****************************************************************************80
%
%% C8 evaluates C function #8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 June 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the evaluation point.
%
%    Output, real VALUE, the value of C(X).
%
  value = 0.0;

  return
end
function value = exact8 ( x )

%*****************************************************************************80
%
%% EXACT8 evaluates exact solution #8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 June 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the evaluation point.
%
%    Output, real VALUE, the value of U(X).
%
  value = x .* ( 1.0 - x ) .* exp ( x )         .* ( x <= 2.0 / 3.0 ) ...
        + x .* ( 1.0 - x ) .* exp ( 2.0 / 3.0 ) .* ( 2.0 / 3.0 < x );

  return
end
function value = exact_ux8 ( x )

%*****************************************************************************80
%
%% EXACT_UX8 evaluates the derivative of exact solution #8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 June 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the evaluation point.
%
%    Output, real VALUE, the value of dUdX(X).
%
  value = ( 1.0 - x - x .* x ) .* exp ( x )         .* ( x <= 2.0 / 3.0 ) ...
        + ( 1.0 - 2.0 * x )    .* exp ( 2.0 / 3.0 ) .* (      2.0 / 3.0 < x );

  return
end
function value = f8 ( x )

%*****************************************************************************80
%
%% F8 evaluates right hand side function #8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 June 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the evaluation point.
%
%    Output, real VALUE, the value of F(X).
%
  value = x .* ( x + 3.0 ) .* exp ( x )         .* ( x <= 2.0 / 3.0 ) ...
        +               2.0 * exp ( 2.0 / 3.0 ) .* (      2.0 / 3.0 < x );

  return
end

