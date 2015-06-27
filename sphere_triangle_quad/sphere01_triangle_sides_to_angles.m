function [ a, b, c ] = sphere01_triangle_sides_to_angles ( as, bs, cs )

%*****************************************************************************80
%
%% SPHERE01_TRIANGLE_SIDES_TO_ANGLES computes triangle angles on the unit sphere.
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
%  Parameters:
%
%    Input, real AS, BS, CS, the (geodesic) length of the 
%    sides of the triangle.
%
%    Output, real A, B, C, the spherical angles of the triangle.
%    Angle A is opposite the side of length AS, and so on.
%
  asu = as;
  bsu = bs;
  csu = cs;
  ssu = ( asu + bsu + csu ) / 2.0;

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / ...
                  ( sin ( ssu ) * sin ( ssu - asu )     ) );

  a = 2.0 * atan ( tan_a2 );

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / ...
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) );

  b = 2.0 * atan ( tan_b2 );

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / ...
                  ( sin ( ssu ) * sin ( ssu - csu )     ) );

  c = 2.0 * atan ( tan_c2 );

  return
end
