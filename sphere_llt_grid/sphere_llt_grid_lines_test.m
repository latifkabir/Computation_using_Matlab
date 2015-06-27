function sphere_llt_grid_lines_test ( )

%*****************************************************************************80
%
%% SPHERE_LLT_GRID_LINES_TEST tests SPHERE_LLT_GRID_LINES.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 April 2015
%
%  Author:
%
%    John Burkardt
%
  lat_num = 3;
  long_num = 4;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'SPHERE_LLT_GRID_LINES_TEST\n' );
  fprintf ( 1, '  SPHERE_LLT_GRID_LINES computes grid lines\n' );
  fprintf ( 1, '  on a sphere in 3D.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of latitudes is  %d\n', lat_num );
  fprintf ( 1, '  Number of longitudes is %d\n', long_num );

  line_num = sphere_llt_grid_line_count ( lat_num, long_num );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of line segments is %d\n', line_num );

  line = sphere_llt_grid_lines ( lat_num, long_num, line_num );

  i4mat_print ( line_num, 2, line, '  Grid line vertices:' );

  return
end
