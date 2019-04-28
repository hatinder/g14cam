function [ u, v, p ] = uvp_poiseuille ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% UVP_POISEUILLE returns the Poiseuille solution.
%
%  Location:
%
%    http://people.sc.fsu.edu/~jburkardt/m_src/navier_stokes_2d_exact/uvp_poiseuille.m
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 July 2015
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real NU, the kinematic viscosity.
%
%    Input, real RHO, the fluid density.
%
%    Input, integer N, the number of nodes.
%
%    Input, real X(N), Y(N), the coordinates of nodes.
%
%    Input, real T, the current time.
%
%    Output, real U(N), V(N), the X and Y velocity.
%
%    Output, real P(N), the pressure.
%
  u(1:n) = 1.0 - y(1:n) .^ 2;

  v(1:n) = 0.0;

  p(1:n) = -2.0 * rho * nu * x(1:n);

  return
end

