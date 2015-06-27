function sparse_count_test07 ( dim_min, dim_max, level_max_min, ...
  level_max_max )

%*****************************************************************************80
%
%% SPARSE_COUNT_TEST07 tests ONN_L_SIZE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    15 January 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer DIM_MIN, the minimum spatial dimension to consider.
%
%    Input, integer DIM_MAX, the maximum spatial dimension to consider.
%
%    Input, integer LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
%
%    Input, integer LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SPARSE_COUNT_TEST07\n' );
  fprintf ( 1, '  ONN_L_SIZE returns the number of\n' );
  fprintf ( 1, '  distinct points in an ONN_L sparse grid made from \n' );
  fprintf ( 1, '  product grids formed from open non-nested\n' );
  fprintf ( 1, '  quadrature rules with linear growth, including:\n' );
  fprintf ( 1, '  * LG_L, the Gauss Laguerre Linear Growth Family;\n' );
  fprintf ( 1, '  * GJ_L, the Gauss Jacobi Linear Growth Family;\n' );
  fprintf ( 1, '  * GLG_L, the Generalized Gauss Laguerre Linear Growth Family;\n' );
  fprintf ( 1, '  * GW_L, any Golub Welsch Linear Growth Family;\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '   DIM: ' );

  for dim_num = dim_min : dim_max
    fprintf ( 1, '  %10d', dim_num );
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '   LEVEL_MAX\n' );
  fprintf ( 1, '\n' );

  for level_max = level_max_min : level_max_max
    fprintf ( 1, '    %4d', level_max );
    for dim_num = dim_min : dim_max
      point_num = onn_l_size ( dim_num, level_max );
      fprintf ( 1, '  %10d', point_num );
    end
    fprintf ( 1, '\n' );
  end

  return
end
