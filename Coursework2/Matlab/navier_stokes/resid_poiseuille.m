function [ ur, vr, pr ] = resid_poiseuille ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% RESID_POISEUILLE evaluates the Poiseuille residual.
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
%    Output, real UR(N), VR(N), PR(N), the residuals.
%

%
%  Get the right hand side functions.
%
  [ f, g, h ] = rhs_poiseuille ( nu, rho, n, x, y, t );
%
%  Form the functions and derivatives for the left hand side.
%
  u(1:n) = 1.0 - y(1:n) .^ 2;
  dudt(1:n) = 0.0;
  dudx(1:n) = 0.0;
  dudxx(1:n) = 0.0;
  dudy(1:n) = - 2.0 * y(1:n);
  dudyy(1:n) = -2.0;

  v(1:n) = 0.0;
  dvdt(1:n) = 0.0;
  dvdx(1:n) = 0.0;
  dvdxx(1:n) = 0.0;
  dvdy(1:n) = 0.0;
  dvdyy(1:n) = 0.0;

  p(1:n) = - 2.0 * rho * nu * x(1:n);
  dpdx(1:n) = - 2.0 * rho * nu;
  dpdy(1:n) = 0.0;

  ur(1:n) = dudt(1:n) - nu * ( dudxx(1:n) + dudyy(1:n) ) ...
    + u(1:n) .* dudx(1:n) + v(1:n) .* dudy(1:n) + dpdx(1:n) / rho - f(1:n);

  vr(1:n) = dvdt(1:n) - nu * ( dvdxx(1:n) + dvdyy(1:n) ) ...
    + u(1:n) .* dvdx(1:n) + v(1:n) .* dvdy(1:n) + dpdy(1:n) / rho - g(1:n);

  pr(1:n) = dudx(1:n) + dvdy(1:n) - h(1:n);

  return
end
