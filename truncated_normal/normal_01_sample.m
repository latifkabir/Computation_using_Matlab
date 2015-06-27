function [ x, seed ] = normal_01_sample ( seed )

%*****************************************************************************80
%
%% NORMAL_01_SAMPLE samples the standard normal probability distribution.
%
%  Discussion:
%
%    The standard normal probability distribution function (PDF) has
%    mean 0 and standard deviation 1.
%
%  Method:
%
%    The Box-Muller method is used, which is efficient, but
%    generates two values at a time.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the standard normal PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  [ r1, seed ] = r8_uniform_01 ( seed );
  [ r2, seed ] = r8_uniform_01 ( seed );

  x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );

  return
end
