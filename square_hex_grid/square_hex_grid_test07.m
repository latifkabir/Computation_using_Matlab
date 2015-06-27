function square_hex_grid_test07 ( )

%*****************************************************************************80
%
%% SQUARE_HEX_GRID_TEST07 tests HEX_GRID_LAYERS.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 January 2009
%
%  Author:
%
%    John Burkardt
%
  test_num = 15;

  nodes_per_layer_test = [ ...
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 ];

  fprintf ( 1, '\n' );
  fprintf ( 1, 'SQUARE_HEX_GRID_TEST07\n' );
  fprintf ( 1, '  For a hexagonal grid of points in a coordinate box,\n' );
  fprintf ( 1, '  given NODES_PER_LAYER, the number of grid points\n' );
  fprintf ( 1, '  along the first layer,\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  HEX_GRID_LAYERS computes LAYERS, the number of\n' );
  fprintf ( 1, '  layers.\n' );

  box(1:2,1:2) = [ 1.0, 2.0; ...
                  4.0, 7.0 ];

  box_print_2d ( box );

  fprintf ( 1, '\n' );
  fprintf ( 1, '   NODES  LAYERS\n' );
  fprintf ( 1, '     PER\n' );
  fprintf ( 1, '   LAYER\n' );
  fprintf ( 1, '\n' );

  for test = 1 : test_num
    nodes_per_layer = nodes_per_layer_test ( test );
    layers = hex_grid_layers ( nodes_per_layer, box );
    fprintf ( 1, '  %6d  %6d\n', nodes_per_layer, layers );
  end

  return
end
