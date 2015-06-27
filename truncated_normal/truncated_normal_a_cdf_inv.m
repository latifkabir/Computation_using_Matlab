function x = truncated_normal_a_cdf_inv ( cdf, mu, sigma, a )

%*****************************************************************************80
%
%% TRUNCATED_NORMAL_A_CDF_INV inverts the lower truncated Normal CDF.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 August 2013
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real CDF, the value of the CDF.
%    0.0D+00 <= CDF <= 1.0.
%
%    Input, real MU, SIGMA, the mean and standard deviation of the
%    parent Normal distribution.
%
%    Input, real A, the lower truncation limit.
%
%    Output, real X, the corresponding argument.
%
  if ( cdf < 0.0 || 1.0 < cdf )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRUNCATED_NORMAL_A_CDF_INV - Fatal error!\n' );
    fprintf ( 1, '  CDF < 0 or 1 < CDF.\n' );
    error ( 'TRUNCATED_NORMAL_A_CDF_INV - Fatal error!' )
  end

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  xi_cdf = ( 1.0 - alpha_cdf ) * cdf + alpha_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return
end
