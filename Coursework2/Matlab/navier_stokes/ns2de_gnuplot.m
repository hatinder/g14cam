function ns2de_gnuplot ( header, n, x, y, u, v, p, s )

%*****************************************************************************80
%
%% NS2DE_GNUPLOT writes the Navier-Stokes velocity field to files for GNUPLOT.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    25 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string HEADER, a header to be used to name the files.
%
%    Input, integer N, the number of evaluation points.
%
%    Input, real X(N), Y(N), the coordinates of the evaluation points.
%
%    Input, real U(N), V(N), P(N), the flow field samples.
%
%    Input, real S, a scale factor for the velocity vectors.
%

%
%  Write the data file.
%
  data_filename = strcat ( header, '_data.txt' );

  data_unit = fopen ( data_filename, 'wt' );

  for i = 1 : n
    fprintf ( data_unit, '  %g  %g  %g  %g  %g  %g  %g\n', ...
      x(i), y(i), u(i), v(i), s * u(i), s * v(i), p(i) );
  end

  fclose ( data_unit );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Data written to "%s".\n', data_filename );
%
%  Write the command file.
%
  command_filename = strcat ( header, '_commands.txt' );
  plot_filename = strcat ( header, '.png' );

  command_unit = fopen ( command_filename, 'wt' );

  fprintf ( command_unit, '#  %s\n', command_filename );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, 'set term png\n' );
  fprintf ( command_unit, 'set output "%s"\n', plot_filename );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, '#  Add titles and labels.\n' );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, 'set xlabel "<--- X --->"\n' );
  fprintf ( command_unit, 'set ylabel "<--- Y --->"\n' );
  fprintf ( command_unit, 'set title "Navier-Stokes velocity field"\n' );
  fprintf ( command_unit, 'unset key\n' );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, '#  Add grid lines.\n' );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, 'set grid\n' );
  fprintf ( command_unit, 'set size ratio -1\n' );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, '#  Timestamp the plot.\n' );
  fprintf ( command_unit, '#\n' );
  fprintf ( command_unit, 'set timestamp\n' );
  fprintf ( command_unit, 'plot "%s" using 1:2:5:6 with vectors \\\n', ...
    data_filename );
  fprintf ( command_unit, '  head filled lt 2 linecolor rgb "blue"\n' );
  fprintf ( command_unit, 'quit\n' );

  fclose ( command_unit );

  fprintf ( 1, '  Commands written to "%s".\n', command_filename );

  return
end

