function [ ur, vr, pr ] = resid_vortex ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% RESID_VORTEX evaluates the Vortex residual.
%
%  Location:
%
%    http://people.sc.fsu.edu/~jburkardt/m_src/navier_stokes_2d_exact/resid_vortex.m
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
%    Output, real UR(N), VR(N), PR(N), the residuals in the U, V and P equations.
%

%
%  Get the right hand sides.
%
  [ f, g, h ] = rhs_vortex ( nu, rho, n, x, y, t );
%
%  Make some temporaries.
%
  cx = cos ( pi * x );
  cy = cos ( pi * y );

  sx = sin ( pi * x );
  sy = sin ( pi * y );

  c2x = cos ( 2.0 * pi * x );
  c2y = cos ( 2.0 * pi * y );

  s2x = sin ( 2.0 * pi * x );
  s2y = sin ( 2.0 * pi * y );
%
%  Form the functions and derivatives for the left hand side.
%
  u    =  -                      cx .* sy;
  dudx =                    pi * sx .* sy;
  dudxx =              pi * pi * cx .* sy;
  dudy =  -                 pi * cx .* cy;
  dudyy =              pi * pi * cx .* sy;
  dudt =  zeros ( size ( x ) );

  v    =                         sx .* cy;
  dvdx =                    pi * cx .* cy;
  dvdxx = -            pi * pi * sx .* cy;
  dvdy =  -                 pi * sx .* sy;
  dvdyy = -            pi * pi * sx .* cy;
  dvdt =  zeros ( size ( x ) );

  p =     - 0.25 * rho *    ( c2x + c2y );
  dpdx =  + 0.5  * rho * pi * s2x;
  dpdy =  + 0.5  * rho * pi       * s2y;
%
%  Evaluate the residuals.
%
  ur = dudt + u .* dudx + v .* dudy + ( 1.0 / rho ) * dpdx ...
     - nu * ( dudxx + dudyy ) - f;

  vr = dvdt + u .* dvdx + v .* dvdy + ( 1.0 / rho ) * dpdy ...
     - nu * ( dvdxx + dvdyy ) - g;

  pr = dudx + dvdy - h;

  return
end

