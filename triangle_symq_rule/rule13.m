function [ x, y, w ] = rule13 ( )

%*****************************************************************************80
%
%% RULE13 returns the rule of degree 13.
%
%  Discussion:
%
%    Order 13 (37 pts)
%    1/6 data for 13-th order quadrature with 10 nodes.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 June 2014
%
%  Author:
%
%    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
%    This MATLAB version by John Burkardt.
%
%  Parameters:
%
%    Output, real X(*), Y(*), the coordinates of the nodes.
%
%    Output, real W(*), the weights.
%
  x = [ ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
     -0.39685461846296817096661395897363, ...
      0.00000000000000000000000000000000, ...
     -0.36877346231328110106712038942974, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
     -0.72443410422579609569939260421083, ...
      0.00000000000000000000000000000000 ];
  y = [ ... ...
     -0.56396461592102123502624022053391, ...
     -0.47207168193213434618577598142010, ...
      0.35411102701218903723318259570554, ...
     -0.54446208086261457427222068671521, ...
      0.00000000000000000000000000000000, ...
     -0.40806649765315498179968814806897, ...
     -0.28109188226279360944191830212271, ...
      0.76131746182322280252086625652347, ...
     -0.53930344496737094791136532996013, ...
      0.10684585018901699613124809207464E+01 ];
  w = [ ... ...
      0.65418593445945714253573279005246E-02, ...
      0.21571270093488444532084979604155E-01, ...
      0.30310770119514528295026716361971E-01, ...
      0.23854497740070562467907639645715E-01, ...
      0.11323203959116968208057898102191E-01, ...
      0.48973694128817658616454407740346E-01, ...
      0.30892926213314122881388117756484E-01, ...
      0.20335382082811117457514058128897E-01, ...
      0.20258422938614600267531787728451E-01, ...
      0.52836422050728359852135506640993E-02 ];

  return
end
