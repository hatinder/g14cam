function [ f, g, h ] = rhs_vortex ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% RHS_VORTEX evaluates the Vortex right hand side.
%
%  Location:
%
%    http://people.sc.fsu.edu/~jburkardt/m_src/navier_stokes_2d_exact/rhs_vortex.m
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real NU, the kinematic viscosity.
%
%    Input, real RHO, the density.
%
%    Input, integer N, the number of points at which the solution is to
%    be evaluated.
%
%    Input, real X(N), Y(N), the coordinates of the points.
%
%    Input, real T or T(N), the time coordinate or coordinates.
%
%    Output, real F(N), G(N), H(N), the right hand sides in the U, V and P 
%    equations.
%
  [ r, c ] = size ( x );

  f = - 2.0 * nu * pi ^ 2 * ( cos ( pi * x ) .* sin ( pi * y ) );
  g =   2.0 * nu * pi ^ 2 * ( sin ( pi * x ) .* cos ( pi * y ) );
  h = zeros ( r, c );

  return
end

