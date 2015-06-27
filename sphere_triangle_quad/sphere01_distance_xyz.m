function dist = sphere01_distance_xyz ( xyz1, xyz2 )

%*****************************************************************************80
%
%% SPHERE01_DISTANCE_XYZ computes great circle distances on a sphere.
%
%  Discussion:
%
%    XYZ coordinates are used.
%
%    We assume the points XYZ1 and XYZ2 lie on the same sphere.
%
%    This computation is a special form of the Vincenty formula.
%    It should be less sensitive to errors associated with very small
%    or very large angular separations.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 September 2010
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    "Great-circle distance",
%    Wikipedia.
%
%  Parameters:
%
%    Input, real XYZ1(3), the coordinates of the first point.
%
%    Input, real XYZ2(3), the coordinates of the second point.
%
%    Output, real DIST, the great circle distance between the points.
%
  lat1 = r8_asin ( xyz1(3) );
  lon1 = r8_atan ( xyz1(2), xyz1(1) );

  lat2 = r8_asin ( xyz2(3) );
  lon2 = r8_atan ( xyz2(2), xyz2(1) );

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) ).^2 ...
      + ( cos ( lat1 ) * sin ( lat2 ) ...
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) ).^2;

  top = sqrt ( top );

  bot = sin ( lat1 ) * sin ( lat2 ) ...
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 );

  dist = atan2 ( top, bot );

  return
end
