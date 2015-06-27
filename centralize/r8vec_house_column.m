function v = r8vec_house_column ( n, a, k )

%*****************************************************************************80
%
%% R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
%
%  Discussion:
%
%    The routine returns a vector V that defines a Householder
%    premultiplier matrix H(V) that zeros out the subdiagonal entries of
%    column K of the matrix A.
%
%       H(V) = I - 2 * v * v'
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 April 2013
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the order of the matrix A.
%
%    Input, real A(N,1), column K of the matrix A.
%
%    Input, integer K, the column of the matrix to be modified.
%
%    Output, real V(N,1), a vector of unit L2 norm which defines an
%    orthogonal Householder premultiplier matrix H with the property
%    that the K-th column of H*A is zero below the diagonal.
%
  v(1:n,1) = 0.0;

  if ( k < 1 || n <= k )
    return
  end

  s = sqrt ( sum ( a(k:n,1).^2 ) );

  if ( s == 0.0 )
    return
  end

  if ( a(k,1) < 0.0 )
    v(k,1) = a(k,1) - abs ( s ) ;
  else
    v(k,1) = a(k,1) + abs ( s );
  end
  v(k+1:n,1) = a(k+1:n,1);

  v(k:n,1) = v(k:n,1) / sqrt ( sum ( v(k:n,1).^2 ) );

  return
end
