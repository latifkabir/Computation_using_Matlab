function p4 = angle_half ( p1, p2, p3 )

%*****************************************************************************80
%
%% ANGLE_HALF finds half an angle.
%
%  Discussion:
%
%    The original angle is defined by the sequence of points P1, P2 and P3.
%
%    The point P4 is calculated so that:
%
%      (P1,P2,P4) = (P1,P2,P3) / 2
%
%        P1
%        /
%       /   P4
%      /  .
%     / .
%    P2--------->P3
%
%    Thanks to Cesar Fraga Bobis for pointing out a typographical error in
%    a previous version of this routine.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 December 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real P1(2,1), P2(2,1), P3(2,1), points defining the angle.
%
%    Input, real P4(2,1), a point defining the half angle.
%    The vector P4 - P2 will have unit norm.
%
  p4(1:2,1) = 0.5 * ( ...
      ( p1(1:2,1) - p2(1:2,1) ) / sqrt ( sum ( ( p1(1:2,1) - p2(1:2,1) ).^2 ) ) ...
    + ( p3(1:2,1) - p2(1:2,1) ) / sqrt ( sum ( ( p3(1:2,1) - p2(1:2,1) ).^2 ) ) );

  p4(1:2,1) = p2(1:2,1) + p4(1:2,1) / sqrt ( sum ( ( p4(1:2,1) ).^2 ) );

  return
end
