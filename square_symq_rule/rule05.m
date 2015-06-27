function [ x, w ] = rule05 ( n )

%*****************************************************************************80
%
%% RULE05 returns the rule of degree 5.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    01 July 2014
%
%  Author:
%
%    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
%    This MATLAB version by John Burkardt.
%
%  Reference:
%
%    Hong Xiao, Zydrunas Gimbutas,
%    A numerical algorithm for the construction of efficient quadrature
%    rules in two and higher dimensions,
%    Computers and Mathematics with Applications,
%    Volume 59, 2010, pages 663-676.
%
%  Parameters:
%
%    Input, integer N, the number of nodes.
%
%    Output, real X(2,N), the coordinates of the nodes.
%
%    Output, real W(N), the weights.
%
  xs = [ ...
    0.1775868202077551E-01, ...
   -0.1775868202077539E-01, ...
    0.7788710544649639E+00, ...
   -0.7788710544649639E+00, ...
   -0.7703781288541645E+00, ...
    0.7703781288541645E+00, ...
   -0.7490353914168658E-33 ];
  ys = [ ...
   -0.9659285494001192E+00, ...
    0.9659285494001192E+00, ...
   -0.5715708301251639E+00, ...
    0.5715708301251639E+00, ...
   -0.5829672991828014E+00, ...
    0.5829672991828014E+00, ...
    0.1356144833394667E-33 ];
  ws = [ ...
    0.2246199725165690E+00, ...
    0.2246199725165690E+00, ...
    0.3901817339168917E+00, ...
    0.3901817339168917E+00, ...
    0.3953508381187504E+00, ...
    0.3953508381187504E+00, ...
    0.8081220356417684E+00 ];

  x(1,1:n) = xs(1:n);
  x(2,1:n) = ys(1:n);
  w(1:n) = ws(1:n);

  return
end