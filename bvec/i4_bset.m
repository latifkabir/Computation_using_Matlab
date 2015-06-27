function value = i4_bset ( i4, pos )

%*****************************************************************************80
%
%% I4_BSET returns a copy of an I4 in which the POS-th bit is set to 1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Military Standard 1753,
%    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
%    9 November 1978.
%
%  Parameters:
%
%    Input, integer I4, the integer to be tested.
%
%    Input, integer POS, the bit position, between 0 and 31.
%
%    Output, integer VALUE, a copy of I4, but with the POS-th bit
%    set to 1.
%
  i4_huge = 2147483647;

  value = round ( i4 );

  if ( pos < 0 )

  elseif ( pos < 31 )

    add = 1;

    if ( 0 <= i4 )
      j = i4;
    else
      j = ( i4_huge + i4 ) + 1;
    end

    for k = 1 : pos
      j = floor ( j / 2 );
      add = add * 2;
    end

    if ( mod ( j, 2 ) == 0 )
      value = i4 + add;
    end

  elseif ( pos == 31 )

    if ( 0 < i4 )
      value = - ( i4_huge - i4 ) - 1;
    end

  elseif ( 31 < pos )

    value = i4;

  end

  return
end
