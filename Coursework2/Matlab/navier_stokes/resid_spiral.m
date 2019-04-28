function [ ur, vr, pr ] = resid_spiral ( nu, rho, n, x, y, t )

%*****************************************************************************80
%
%% RESID_SPIRAL evaluates the residuals of the spiral flow problem.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 January 2015
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Maxim Olshanskii, Leo Rebholz,
%    Application of barycenter refined meshes in linear elasticity
%    and incompressible fluid dynamics,
%    ETNA: Electronic Transactions in Numerical Analysis, 
%    Volume 38, pages 258-274, 2011.
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
  [ f, g, h ] = rhs_spiral ( nu, rho, n, x, y, t );
%
%  Form the functions and derivatives for the left hand side.
%
  u(1:n) = ( 1.0 + nu * t ) * 2.0 ...
    .* ( x(1:n) .^ 4 - 2.0 * x(1:n) .^ 3 + x(1:n) .^ 2 ) ...
    .* ( 2.0 * y(1:n) .^ 3 - 3.0 * y(1:n) .^ 2 + y(1:n) );

  dudt(1:n) = nu * 2.0 ...
    .* ( x(1:n) .^ 4 - 2.0 * x(1:n) .^ 3 + x(1:n) .^ 2 ) ...
    .* ( 2.0 * y(1:n) .^ 3  - 3.0 * y(1:n) .^ 2 + y(1:n) );

  dudx(1:n) = ( 1.0 + nu * t ) * 2.0 ...
    .* ( 4.0 * x(1:n) .^ 3 - 6.0 * x(1:n) .^ 2 + 2.0 * x(1:n) ) ...
    .* ( 2.0 * y(1:n) .^ 3 - 3.0 * y(1:n) .^ 2 + y(1:n) );

  dudxx(1:n) = ( 1.0 + nu * t ) * 2.0 ...
    .* ( 12.0 * x(1:n) .^ 2 - 12.0 * x(1:n) + 2.0 ) ...
    .* ( 2.0 * y(1:n) .^ 3 - 3.0 * y(1:n) .^ 2 + y(1:n) );

  dudy(1:n) = ( 1.0 + nu * t ) * 2.0 ...
    .* ( x(1:n) .^ 4 - 2.0 * x(1:n) .^ 3 + x(1:n) .^ 2 ) ...
    .* ( 6.0 * y(1:n) .^ 2  - 6.0 * y(1:n) + 1.0 );

  dudyy(1:n) = ( 1.0 + nu * t ) * 2.0 ...
    .* ( x(1:n) .^ 4 - 2.0 * x(1:n) .^ 3 + x(1:n) .^ 2 ) ...
    .* ( 12.0 * y(1:n) - 6.0 );

  v(1:n) = - ( 1.0 + nu * t ) * 2.0 ...
    .* ( 2.0 * x(1:n) .^ 3 - 3.0 * x(1:n) .^ 2 + x(1:n) ) ...
    .* ( y(1:n) .^ 4 - 2.0 * y(1:n) .^ 3 + y(1:n) .^ 2 );

  dvdt(1:n) = - nu * 2.0 ...
    .* ( 2.0 * x(1:n) .^ 3 - 3.0 * x(1:n) .^ 2 + x(1:n) ) ...
    .* ( y(1:n) .^ 4 - 2.0 * y(1:n) .^ 3 + y(1:n) .^ 2 );

  dvdx(1:n) = - ( 1.0 + nu * t ) * 2.0 ...
    .* ( 6.0 * x(1:n) .^ 2 - 6.0 * x(1:n) + 1.0 ) ...
    .* ( y(1:n) .^ 4 - 2.0 * y(1:n) .^ 3 + y(1:n) .^ 2 );

  dvdxx(1:n) = - ( 1.0 + nu * t ) * 2.0 ...
    .* ( 12.0 * x(1:n) - 6.0 ) ...
    .* ( y(1:n) .^ 4 - 2.0 * y(1:n) .^ 3 + y(1:n) .^ 2 );

  dvdy(1:n) = - ( 1.0 + nu * t ) * 2.0 ...
    .* ( 2.0 * x(1:n) .^ 3 - 3.0 * x(1:n) .^ 2 + x(1:n) ) ...
    .* ( 4.0 * y(1:n) .^ 3 - 6.0 * y(1:n) .^ 2 + 2.0 * y(1:n) );

  dvdyy(1:n) = - ( 1.0 + nu * t ) * 2.0 ...
    .* ( 2.0 * x(1:n) .^ 3 - 3.0 * x(1:n) .^ 2 + x(1:n) ) ...
    .* ( 12.0 * y(1:n) .^ 2 - 12.0 * y(1:n) + 2.0 );

  p(1:n) = rho * y(1:n);
  dpdx(1:n) = 0.0;
  dpdy(1:n) = rho;

  ur(1:n) = dudt(1:n) - nu * ( dudxx(1:n) + dudyy(1:n) ) ...
    + u(1:n) .* dudx(1:n) + v(1:n) .* dudy(1:n) + dpdx(1:n) / rho - f(1:n);

  vr(1:n) = dvdt(1:n) - nu * ( dvdxx(1:n) + dvdyy(1:n) ) ...
    + u(1:n) .* dvdx(1:n) + v(1:n) .* dvdy(1:n) + dpdy(1:n) / rho - g(1:n);

  pr(1:n) = dudx(1:n) + dvdy(1:n) - h(1:n);

  return
end
