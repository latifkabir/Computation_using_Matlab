function value = f4 ( x )

%*****************************************************************************80
%
%% F4 is an example of a function whose "area under the curve" is to be plotted.
%
%  Discussion:
%
%    This function should be written in such a way that it can accept
%    an input X that is a vector.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 February 2013
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument.
%
%    Output, real VALUE, the value of the function.
%
  value = exp ( - 0.5 * ( x - 5 ).^2 ) / sqrt ( 2.0 * pi );

  return
end

