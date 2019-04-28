function [ f, g, h ] = rhs_poiseuille ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% RHS_POISEUILLE evaluates the Poiseuille right hand side.
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
%    Output, real F(N), G(N), H(N), the right hand sides.
%
  f(1:n) = 0.0;

  g(1:n) = 0.0;

  h(1:n) = 0.0;

  return
end
