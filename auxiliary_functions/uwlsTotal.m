function [ u, v, dopx, dopy, hdop ] = uwlsTotal( rSpeed, rHeading )
% UWLSTOTAL Computes a total velocity from radials using least-squares
%
%    [ u, v, dopx, dopy, hdop ] = uwlsTotal( rSpeed, rHeading ) returns the
%    eastward (u) and northward (v) components of the total velocity
%    resulting from the unweighted least squares fit of radial velocity
%    measurements provided by rSpeed and rHeading.  The horizontal dilution
%    of precision (hdop) and corresponding eastward (dopx) and northward
%    (dopy) components are also returned as a metric for assessing the
%    quality of the solution.
%
%    INPUT:
%        rSpeed   - Column vector (n* x 1) of radial velocity magnitude
%        rHeading - Column vector (n* x 1) of radial velocity heading in
%                   degrees counterclockwise from +x (east)
%
%        *n >= 2
%
%    OUTPUT:
%        u    - Eastward total velocity
%        v    - Northward total velocity
%        dopx - Eastward dilution of precision
%        dopy - Northward dilution of precision
%        hdop - Horizontal dilution of precision
%
%    All output values are scalars and the units of eastward and northward
%    velocity are the same as the units for rSpeed.
%
%    Radial velocities are related to the total velocity by projection of
%    the eastward and northward components of the total velocity onto the
%    radial heading:
%
%        rSpeed = u*cosd(rHeading) + v*sind(rHeading)
%
%    and in matrix form as
%
%        rSpeed = X*b
%
%    where
%
%        X = [ cosd(rHeading) sind(rHeading) ]
% b = [ u; v ] %
%    The unweighted least squares solution for b is given by
%
%        b = inv( X'*X ) * X' * rSpeed
%
%    The geometric portion of the covariance matrix for b is
%
%        C = inv( X'*X )
%
%    which provides the dilution of precision along the diagonal where
%
%        dopx = sqrt( C(1,1) )
%        dopy = sqrt( C(2,2) )
%        hdop = sqrt( C(1,1) + C(2,2) )
%
%    References
%
%        Least-squares methods for the extraction of surface currents from
%        CODAR crossed-loop data: Application at ARSLOE, Lipa, B. J., D. E.
%        Barrick, IEEE J. Ocean. Eng., 1983.
%
%        Shipborne measurement of surface current fields by HF radar,
%        Gurgel, K. -W., L'Onde Electr., 1994.
%
%        On the accuracy of HF radar suface current measurements:
%        Intercomparisons with ship-based sensors, Chapman, R. D., L. K.
%        Shay, H. C. Graber, J. B. Edson, A. Karachintsev, C. L. Trump, D.
%        B. Ross, J. Geophys. Res., 1997.
%
%        Guide to GPS Positioning, Wells, D. E., et al., Can. GPS Assoc.,
%        Fredericton, N. B., Canada, 1986.
%
%    Acknowledgement
%
%         This function was adapted from the TOTCALC function (ver. 3.0, Apr
%         2004) included in the HFRadarmap toolbox (ver. 4.1, Sept 2004) written
%         by Mike Cook while at the Naval Postgraduate School Department of
%         Oceanography.
% 
%     HF-Radar Network
%     Scripps Institution of Oceanography
%     Coastal Observing Research and Development Center

X = zeros( numel( rHeading ) , 2 );
X(:,1) = cosd( rHeading );
X(:,2) = sind( rHeading );
C = inv( X' * X );
dopx = sqrt( C(1,1) );
dopy = sqrt( C(2,2) );
hdop = sqrt( C(1,1) + C(2,2) );
b = C * X' * rSpeed;
u = b(1);
v = b(2);
end
% Shared with Luke Kachelein by Mark Otero, May 3, 2021
% [ u, v, dopx, dopy, hdop ] = uwlsTotal([1 1]',[45 135]') % simple example
% that passes a sanity check