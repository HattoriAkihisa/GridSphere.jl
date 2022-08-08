"""
Splits each given triangle into numNewTrianglesPerTriangle equal triangles.
Assumes triangles is a 3-D array. Each row represents a triangle. Columns
represent x, y, and z coordinates. Depth represents the A, B, and C triangle
vertices. Ensures that every vertex is an A, B, and C vertex in at least one
triangle. Work in Cartesian coordinates, projecting onto a unit sphere at each
iteration. The resulting points are evenly spaced in great circle angle, not
latitude and/or longitude, which ensures correct handling of the poles.
numDimensions is the number of dimensions of space. Assumes triangles have
numVerticesInATriangle vertices each. Does not handle arbitrary input for the
last three parameters. Only handles numDimensions = 3,
numVerticesInATriangle = 3, and numNewTrianglesPerTriangle = 4.
-----------------------------------------------------------------------
N. Teanby   13-01-04    Original IDL (Interactive Data Language) code
(available at http://www.atm.ox.ac.uk/user/teanby/software.html#icos)
Kurt von Laven    18-09-13    Ported simplified version to GNU Octave.
-----------------------------------------------------------------------
   usage: subTriangles = DivideTriangles(triangles, numDimensions, ...
                           numVerticesInATriangle, numNewTrianglesPerTriangle)
"""
function dividetriangles(triangles, num_dimensions, num_vertices_in_triangle, num_new_triangles_per_triangle)
    # Get the original triangle vertices
    oldAs = triangles[:, :, 1];
    oldBs = triangles[:, :, 2];
    oldCs = triangles[:, :, 3];

    # Find the midpoints of each side.
    Ps = (oldAs + oldBs) ./ 2;
    Qs = (oldBs + oldCs) ./ 2;
    Rs = (oldCs + oldAs) ./ 2;

    # Normalize midpoints onto the surface of a unit sphere.
    scaling = 1 / norm(Ps[1,:]);
    unitPs = scaling * Ps;
    unitQs = scaling * Qs;
    unitRs = scaling * Rs;

    # Find the sub-triangles' vertices. Ensure that every point gets used as an
    # A, B, and C in at least one triangle.
    newAs = [oldAs; unitPs; unitRs; unitQs];
    newBs = [unitPs; oldBs; unitQs; unitRs];
    newCs = [unitRs; unitQs; oldCs; unitPs];

    subtriangles = reshape(
        [newAs newBs newCs], 
        num_new_triangles_per_triangle * size(triangles, 1),
        num_dimensions,
        num_vertices_in_triangle
    );
end