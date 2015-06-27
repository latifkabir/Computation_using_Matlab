function [ x, w ] = rule08 ( n )

%*****************************************************************************80
%
%% RULE08 returns the rule of degree 8.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 July 2014
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
    -.2272218649369121,0.2786782798124801, ...
    0.9215721988395638,-.5229427015551803, ...
    0.8309170589376613,-.6080254018462903, ...
    -.9822549066167084,0.4959470731361600E-01, ...
    0.5910013957537859,0.3626589212754838, ...
    -.9369162594801185,-.8850131220663160, ...
    -.1934658240272289,0.5772453681919104, ...
    0.9213070164035271,-.7176037958340967 ];
  ys = [ ...
    0.8703146041404044,0.9856262640199153, ...
    0.2224095500358621,-.9282264259882677, ...
    0.8435111761265234,0.5825946042673711, ...
    -.8211266831948021,-.6917239446781449, ...
    -.2614406969784849,0.5198121135620160, ...
    0.2153771996329335,0.9090384216207131, ...
    0.3526321874643216E-01,-.9622555555961493, ...
    -.7082682674817122,-.4130619139730907 ];
  ws = [ ...
    0.1444235515837947,0.5206905878850610E-01, ...
    0.1365819925705312,0.1136963049256808, ...
    0.1156201396846171,0.2194056396025883, ...
    0.4570142629159132E-01,0.3040158377300561, ...
    0.3227772111095287,0.3341175763908440, ...
    0.1202823186503543,0.6155232134515501E-01, ...
    0.4037250536437860,0.8510021531985533E-01, ...
    0.1026971066172272,0.2666613704920739 ];

  x(1,1:n) = xs(1:n);
  x(2,1:n) = ys(1:n);
  w(1:n) = ws(1:n);

  return
end