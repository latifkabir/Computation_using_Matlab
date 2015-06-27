function s_new = s_blank_delete ( s )

%*****************************************************************************80
%
%% S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
%
%  Discussion:
%
%    All TAB characters are also removed.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 March 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character S(*), the string to be transformed.
%
%    Output, character S_NEW(*), the transformed string.
%
  TAB = 9;

  s_len = length ( s );
  s_new = '';

  put = 0;

  for get = 1 : s_len

    if ( s(get) ~= ' ' & s(get) ~= TAB )
      put = put + 1;
      s_new(put) = s(get);
    end

  end

  return
end
