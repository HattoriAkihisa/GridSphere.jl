"""
Converts Cartesian coordinates to latitudes and longitudes.
Cartesian coordinates should lie on a unit sphere. The origin of the
coordinate system is at the center of the sphere. The positive x-axis passes
through 0 degrees longitude, the positive y-axis passes through +90 degrees
longitude, and the positive z-axis passes through +90 degrees latitude.
Latitude and longitude are given in degrees with latitudes of [-90, 90] and
longitudes of (-180, 180]. The point (0, 0, 0) at the center of the sphere is
defined to have a latitude and longitude of NaN (not a number). Refer to
http://www.geomidpoint.com/calculation.html for more information.

   usage: [latsInDegrees, longsInDegrees] = Cartesian2Spherical(Xs, Ys, Zs)
"""
function cartesian2spherical(xs, ys, zs)
    lats = @__dot__ asind(zs ./ sqrt(xs^2 + ys^2 + zs^2));
    
    # Create a mask for points for which division by zero occurred. These
    # points correspond to the center of the sphere.
    # Overwrite whatever bad value indicator (i.e., NaN + NaNi) resided in these
    # points with the more appropriate indiciator, NaN (not a number).
    centers = map(!isfinite, lats)
    lats[centers] .= NaN

    lons = @__dot__ atand(ys, xs)

    return lats, lons
end