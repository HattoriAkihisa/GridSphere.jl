module GridSphere

export gridsphere

using LinearAlgebra

include("divide_triangles.jl")
include("cartesian2spherical.jl")

const golden_ratio = (1 + sqrt(5)) / 2;


"""
Computes the element-wise logarithms of the given values.
bases are the bases for which the logarithms are taken.
   usage: logarithms = Logarithm(values, bases)
"""
logarithm(values, bases) = @__dot__ log(values) / log(bases)

"""
Generates a geodesic grid, a set of nearly evenly spaced points on a sphere.

Points are sorted in order of increasing latitude. Ties between latitudes that
differ only by floating point errors are broken by a secondary sort in order
of increasing longitude. Tesselates a sphere with triangles. Starts with an
icosahedron circumscribed by a unit sphere. Divides the triangles into 4 new
equal triangles. Projects the new vertices onto the unit sphere. Repeats until
reaching the desired resolution. Refer to
http://www.scidacreview.org/0904/images/hardware04.jpg for a graphical
depiction of the algorithm.

-----------------------------------------------------------------------
N. Teanby   13-01-04    Original IDL (Interactive Data Language) code
(available at http://www.atm.ox.ac.uk/user/teanby/software.html#icos)
Kurt von Laven    18-09-13    Ported simplified version to GNU Octave.

-----------------------------------------------------------------------

   usage: latlons = gridsphere(n)

   n = 12, 2 + 10*4^k
"""
function gridsphere(n::Int)
    # Set the golden ratio.
    ϕ = golden_ratio

    # The number of faces in an icosahedron
    num_triangles_in_icosahedron = 20

    # The dimensionality of the sphere.
    num_dimentions = 3

    # The number of vertices in a triangle.
    num_vertices_in_triangle = 3

    # The number of triangle produced
    num_new_triangles_per_ptriangle = 4

    # Find the 12 vertices of an icosahedron inscribed in a unit sphere.
    vertices = [
        0 ϕ 1;
        0 -ϕ 1;
        0 ϕ -1;
        0 -ϕ -1;
        1 0 ϕ;
        -1 0 ϕ;
        1 0 -ϕ;
        -1 0 -ϕ;
        ϕ 1 0;
        -ϕ 1 0;
        ϕ -1 0;
        -ϕ -1 0;
    ] ./ norm([1, ϕ])

    # Find the vertices of the starting triangles. Each row represents a
    # triangle. Columns represent x, y, and z coordinates. Depth represents the
    # A, B, and C triangle vertices. Ensure that each vertex is an A, B, and C
    # vertex in at least one triangle.
    As = [2, 5, 9, 7, 11, 4, 6, 2, 1, 5, 3, 9, 8, 7, 12, 12, 6, 1, 3, 10]
    Bs = [4, 11, 5, 9, 7, 2, 12, 5, 6, 9, 1, 7, 3, 4, 8, 6, 1, 3, 10, 8]
    Cs = [11, 2, 11, 11, 4, 12, 2, 6, 5, 1, 9, 3, 7, 8, 4, 10, 10, 10, 8, 12]

    triangles = reshape(
        [vertices[As,:] vertices[Bs,:] vertices[Cs,:]], 
        num_triangles_in_icosahedron, num_dimentions, num_vertices_in_triangle
    )


    # Determines the number of divisions that will result in a grid with as
    # close to the number of points as the user requested as possible.
    num_divisions = ceil(logarithm((n - 2) / 25, num_new_triangles_per_ptriangle))

    for i in 1:num_divisions
        # Split each triangle into 4 smaller triangles.
        triangles = dividetriangles(
            triangles,
            num_dimentions,
            num_vertices_in_triangle,
            num_new_triangles_per_ptriangle)
    end

    lats, lons = cartesian2spherical(
        triangles[:,1,1],
        triangles[:,2,1],
        triangles[:,3,1])

    latlon = sortslices([lats lons], dims = 1, by = x -> x[1]) |>
        M -> unique(M, dims = 1)

    return latlon
end

end # module
