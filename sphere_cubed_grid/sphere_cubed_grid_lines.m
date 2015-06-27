function line_data = sphere_cubed_grid_lines ( n, line_num )

%*****************************************************************************80
%
%% SPHERE_CUBED_GRID_LINES computes the lines on a cubed sphere grid.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 May 2015
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of sections into which each face of
%    the cube is to be divided.
%
%    Input, integer LINE_NUM, the number of lines.
%
%    Output, real LINE_DATA(LINE_NUM,3,2), for each line I, the X/Y/Z 
%    coordinates of the start and end of a line segment on the grid.
%
  line_data = zeros ( line_num, 3, 2 );

  l = 0;
%
%  If N = 1, the corners form 12 lines.
%
  if ( n == 1 )
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n );
    return
%
%  If 1 < N, each of 8 corners connects to three neighboring edges.
%
  else
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 1, 0, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n-1, 0, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n-1, n, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n-1, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, 1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 1, n, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n-1, 0 );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0 );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 1, 0, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n-1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n-1, 0, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 1, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n-1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n-1, n, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n-1, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, n-1 );

    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 1, n, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n-1, n );
    l = l + 1;
    line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n );
    line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n-1 );
  end
%
%  If 2 < N, then each of the 12 edges includes lines.
%
  if ( 2 < n )

    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i,   0, 0 );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, 0 );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n,   i, 0 );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, 0 );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n-i,   n, 0 );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n-i-1, n, 0 );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i,   0 );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i-1, 0 );
    end

    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i,   0, n );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, n );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n,   i, n );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, n );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n-i,   n, n );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n-i-1, n, n );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i,   n );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i-1, n );
    end

    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, i   );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, i+1 );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, i   );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, 0, i+1 );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, i   );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, n, i+1 );
    end
    for i = 1 : n - 2
      l = l + 1;
      line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, i   );
      line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, n, i+1 );
    end

  end
%
%  Lines that belong to one of the six faces.
%
  if ( 1 < n )
%% 000 : nn0
    for i = 1 : n - 1
      for j = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i, j,   0 );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i, j+1, 0 );
      end
    end
    for j = 1 : n - 1
      for i = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i,   j, 0 );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i+1, j, 0 );
      end
    end
%
%  00n : nnn
%
    for i = 1 : n - 1
      for j = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i, j,   n );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i, j+1, n );
      end
    end
    for j = 1 : n - 1
      for i = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i,   j, n );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i+1, j, n );
      end
    end
%
%  000:n0n
%
    for i = 1 : n - 1
      for j = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i, 0, j   );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i, 0, j+1 );
      end
    end
    for j = 1 : n - 1
      for i = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i,   0, j );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, j );
      end
    end
%
%  0n0:nnn
%
    for i = 1 : n - 1
      for j = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i, n, j   );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i, n, j+1 );
      end
    end
    for j = 1 : n - 1
      for i = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, i,   n, j );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, i+1, n, j );
      end
    end
%
%  000:0nn
%
    for i = 1 : n - 1
      for j = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, i, j   );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, i, j+1 );
      end
    end
    for j = 1 : n - 1
      for i = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, 0, i,   j );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, 0, i+1, j );
      end
    end
%
%  n00:nnn
%
    for i = 1 : n - 1
      for j = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, i, j   );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, i, j+1 );
      end
    end
    for j = 1 : n - 1
      for i = 0 : n - 1
        l = l + 1;
        line_data(l,1:3,1) = sphere_cubed_grid_ijk_to_xyz ( n, n, i,   j );
        line_data(l,1:3,2) = sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, j );
      end
    end

  end

  if ( l ~= line_num )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SPHERE_CUBED_GRID_LINES - Fatal error!\n' );
    fprintf ( 1, '  LINE_NUM = %d\n', line_num );
    fprintf ( 1, '  L = %d\n', l );
    error ( 'SPHERE_CUBED_GRID_LINES - Fatal error!' );
  end

  return
end
