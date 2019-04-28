function [ u, v, p ] = uvp_vortex ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% UVP_VORTEX evaluates the Vortex solution.
%
%  Location:
%
%    http://people.sc.fsu.edu/~jburkardt/m_src/navier_stokes_2d_exact/uvp_vortex.m
%
%  Discussion:
%
%    The given velocity and pressure fields are exact solutions for the 2D 
%    incompressible time-dependent Navier Stokes equations over any region.
%
%    To define a typical problem, one chooses a bounded spatial region 
%    and a starting time, and then imposes boundary and initial conditions
%    by referencing the exact solution appropriately.
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
%    Output, real U(N), V(N), P(N), the velocity components and
%    pressure at each of the points.
%
  cx = cos ( pi * x );
  cy = cos ( pi * y );
  c2x = cos ( 2.0 * pi * x );
  c2y = cos ( 2.0 * pi * y );
  sx = sin ( pi * x );
  sy = sin ( pi * y );

  u = - cx .* sy;
  v =   sx .* cy;
  p = - 0.25 * rho * ( c2x + c2y );

  return
end

