function area = sphere01_area ( )

%*****************************************************************************80
%
%% SPHERE01_AREA returns the area of the surface of the unit sphere in 3D.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    06 January 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Output, real AREA, the area.
%
  r = 1.0;
  area = 4.0 * pi * r * r;

  return
end
