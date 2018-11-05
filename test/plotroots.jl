using VoronoiDelaunay
# import Cairo, Fontconfig

using GRPF

using Gadfly
import Gadfly: RGB
using DataFrames


push!(LOAD_PATH, "C:/dev/LongwaveModeSolver/src")
using ModeFinder
using Geophysics

function groundf(θ)
    ω = 2π*24e3
    k = ω/ModeFinder.speedoflight*1e3

    # Seawater
    σ = 5.
    ϵᵣ = 88.

    electrons = Constituent(-fundamentalcharge, mₑ,
                            h -> waitprofile(h, 75., 0.3), collisionprofile)

    Xsol = ModeFinder.integratethroughionosphere(θ, ω, k, 90., 10., 50.,
                                                 electrons, 0.5e-4, 0., 0., -1.)
    X = Xsol[end]

    N11, N22, D11, D22 = ModeFinder.fresnelgroundreflection(θ, ω, σ, ϵᵣ)

    return (N11 - X[1,1]*D11)*(N22 - X[2,2]*D22) - X[2,1]*X[1,2]*D11*D22
end

function reflectionf(θ)
    ω = 2π*24e3
    k = ω/ModeFinder.speedoflight*1e3

    # Seawater
    σ = 5.
    ϵᵣ = 88.

    topheight = 90.
    bottomheight = 45.
    referenceheight = 50.
    reflectionheight = 68.

    electrons = Constituent(-fundamentalcharge, mₑ,
                            h -> waitprofile(h, 75., 0.3), collisionprofile)

    # Calculate `X` at `bottomheight`
    Xsol = ModeFinder.integratethroughionosphere(θ, ω, k, topheight, bottomheight, referenceheight,
                                                 electrons, 0.5e-4, 0., 0., -1.)
    X = Xsol[end]

    # Refer `X` from `bottomheight` to `referenceheight`
    ModeFinder.integratethroughfreespace!(X, θ, k, bottomheight, reflectionheight, referenceheight)

    # Calculate ground reflection matrix as `N` and `D` referred to `referenceheight`
    # XXX: Blows up at low θ?
    # XXX: MUST BE DUE TO HANKEL FCN XXX
    N11, N22, D11, D22 = ModeFinder.groundreflection(θ, ω, k, σ, ϵᵣ, reflectionheight, referenceheight)

    return (N11 - X[1,1]*D11)*(N22 - X[2,2]*D22) - X[2,1]*X[1,2]*D11*D22
end

function mygroundf(θ)
    ω = 2π*24e3
    k = ω/ModeFinder.speedoflight*1e3

    # Seawater
    σ = 5.
    ϵᵣ = 88.

    electrons = Constituent(-fundamentalcharge, mₑ,
                            h -> waitprofile(h, 75., 0.3), collisionprofile)

    Xsol = ModeFinder.integratethroughionosphere(θ, ω, k, 90., 10., 50.,
                                                 electrons, 0.5e-4, 0., 0., -1.)
    X = Xsol[end]

    N11, N22, D11, D22 = ModeFinder.fresnelgroundreflection(θ, ω, σ, ϵᵣ)

    return (N11 - X[1,1]*D11)*(N22 - X[2,2]*D22) - X[2,1]*X[1,2]*D11*D22
end

"""
Linearly map function values within domain from `min_coord` to `max_coord`.
"""
function mapfunctionval(z, ra, rb, ia, ib)
    zr = ra*real(z) + rb
    zi = ia*imag(z) + ib
    complex(zr, zi)
end
function mapfunctionval!(z, ra, rb, ia, ib)
    for ii in eachindex(z)
        z[ii] = mapfunctionval(z[ii], ra, rb, ia, ib)
    end
end

function geom2fcn(pt::AbstractPoint2D, ra, rb, ia, ib)
    complex((getx(pt) - rb)/ra, (gety(pt) - ib)/ia)
end
geom2fcn(edge::VoronoiDelaunay.DelaunayEdge, ra, rb, ia, ib) = (geom2fcn(geta(edge), ra, rb, ia, ib), geom2fcn(getb(edge), ra, rb, ia, ib))
geom2fcn(pa, pb, ra, rb, ia, ib) = complex((pa-rb)/ra, (pb-ib)/ia)

function groundmain()
    # Zb, Ze = ModeFinder.boundaries(24e3)
    Zb = complex(30., -12.)
    Ze = complex(90., 0.)
    r = 0.5 # initial mesh step
    tolerance = 0.01

    origcoords = GRPF.rectangulardomain(Zb, Ze, r)

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    # Adding extra eps() to handle float error
    ra = ((max_coord-eps())-(min_coord+eps()))/(rmax-rmin)
    rb = (max_coord-eps()) - ra*rmax

    ia = ((max_coord-eps())-(min_coord+eps()))/(imax-imin)
    ib = (max_coord-eps()) - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(10000)

    @assert minimum(real(origcoords)) >= min_coord
    @assert minimum(imag(origcoords)) >= min_coord
    @assert maximum(real(origcoords)) <= max_coord
    @assert maximum(imag(origcoords)) <= max_coord

    tess, 𝓔, quadrants, phasediffs = GRPF.tesselate!(tess, newnodes, pt -> groundf(geom2fcn(pt, ra, rb, ia, ib)),
                                         e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

    𝐶 = GRPF.contouredges(tess, 𝓔)
    regions = GRPF.evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))
    #
    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = GRPF.rootsandpoles(regions,
                                                                                  quadrants,
                                                                                  e -> geom2fcn(e, ra, rb, ia, ib))

    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    # zroots:
    # 18-element Array{Complex{Float64},1}:
    #   73.77652877697841 - 1.5054687500000004im
    #  49.016893627954786 - 4.341517857142856im
    #   79.76506294964031 - 0.940755208333334im
    #  49.491906474820134 - 1.9164062500000008im
    #  55.780125899280584 - 3.37890625im
    #   67.58430755395683 - 1.6302083333333333im
    #  34.373875899280584 - 2.3157552083333317im
    #   62.29856115107913 - 2.7968750000000013im
    #    84.2974370503597 - 0.37760416666666513im
    #   77.61960431654674 - 1.79453125im
    #    83.0328237410072 - 0.9984375000000005im
    #   61.74955035971224 - 1.7632812500000008im
    #   55.84869604316546 - 1.8470052083333337im
    #   42.48156474820143 - 2.0523437500000004im
    #   41.90422661870504 - 5.75703125im
    #  34.518884892086334 - 7.809375000000001im
    #   68.31497302158273 - 2.652473958333334im
    #   72.98875899280574 - 2.5583333333333327im
    #
    # zpoles:
    # 16-element Array{Complex{Float64},1}:
    #   35.90490107913669 - 9.985677083333334im
    #   72.11330935251797 - 5.096875im
    #    56.6449340527578 - 6.4288194444444455im
    #  54.331160071942435 - 9.238932291666666im
    #   60.25359712230217 - 7.748437500000001im
    #   62.56969424460432 - 5.914843749999999im
    #   78.05193345323741 - 3.6490885416666665im
    #   83.73732831916286 - 1.5759943181818168im
    #   83.92027449468999 - 1.7349950396825389im
    #   50.37005395683454 - 7.178125000000001im
    #   78.61061151079137 - 3.8860677083333335im
    #    68.1946942446043 - 5.57890625im
    #   66.18884892086332 - 6.396875im
    #   43.56789568345322 - 8.274739583333332im
    #  48.419514388489205 - 10.82890625im
    #   73.42176258992805 - 5.1203125im


    edges = [e for e in delaunayedges(tess)]
    df = prepareplotdata(edges, quadrants, phasediffs, (ea, eb) -> geom2fcn(ea, eb, ra, rb, ia, ib))

    p = plot(df, x=:x, y=:y, color=:c, Geom.path, Scale.color_discrete,
             # Coord.cartesian(xmin=60., xmax=70., ymin=-7., ymax=-5.),
             Guide.colorkey(title="𝑓(θ)",
                            labels=["Re > 0, Im ≥ 0", "Re ≤ 0, Im > 0",
                                    "Re &lt; 0, Im ≤ 0", "Re ≥ 0, Im &lt; 0",
                                    "Candidate edges"]),
             Scale.color_discrete_manual(RGB(0.878, 0.749, 0.094), RGB(0.247, 0.125, 0.604),
                                         RGB(0.063, 0.588, 0.404), RGB(0.878, 0.349, 0.094),
                                         RGB(0., 0., 0.), RGB(0.5, 0.5, 0.5)),
             Guide.xlabel("Re(θ) [deg]"), Guide.ylabel("Im(θ) [deg]"))

    draw(SVG("naa_atground.svg", 13inch, 8inch), p)
end


function reflectionmain()
    # Zb, Ze = ModeFinder.boundaries(24e3)
    Zb = complex(40., -9.)
    Ze = complex(90., 0.)
    r = 0.5 # initial mesh step
    tolerance = 0.01

    origcoords = GRPF.rectangulardomain(Zb, Ze, r)

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    # Adding extra eps() to handle float error
    ra = ((max_coord-eps())-(min_coord+eps()))/(rmax-rmin)
    rb = (max_coord-eps()) - ra*rmax

    ia = ((max_coord-eps())-(min_coord+eps()))/(imax-imin)
    ib = (max_coord-eps()) - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(10000)

    @assert minimum(real(origcoords)) >= min_coord
    @assert minimum(imag(origcoords)) >= min_coord
    @assert maximum(real(origcoords)) <= max_coord
    @assert maximum(imag(origcoords)) <= max_coord

    tess, 𝓔, quadrants, phasediffs = GRPF.tesselate!(tess, newnodes, pt -> reflectionf(geom2fcn(pt, ra, rb, ia, ib)),
                                                     e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

    𝐶 = GRPF.contouredges(tess, 𝓔)
    regions = GRPF.evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))
    #
    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = GRPF.rootsandpoles(regions,
                                                                                  quadrants,
                                                                                  e -> geom2fcn(e, ra, rb, ia, ib))

    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    edges = [e for e in delaunayedges(tess)]
    df = prepareplotdata(edges, quadrants, phasediffs, (ea, eb) -> geom2fcn(ea, eb, ra, rb, ia, ib))

    p = plot(df, x=:x, y=:y, color=:c, Geom.path, Scale.color_discrete,
             # Coord.cartesian(xmin=60., xmax=70., ymin=-7., ymax=-5.),
             Guide.colorkey(title="𝑓(θ)",
                            labels=["Re > 0, Im ≥ 0", "Re ≤ 0, Im > 0",
                                    "Re &lt; 0, Im ≤ 0", "Re ≥ 0, Im &lt; 0",
                                    "Candidate edges"]),
             Scale.color_discrete_manual(RGB(0.878, 0.749, 0.094), RGB(0.247, 0.125, 0.604),
                                         RGB(0.063, 0.588, 0.404), RGB(0.878, 0.349, 0.094),
                                         RGB(0., 0., 0.), RGB(0.5, 0.5, 0.5)),
             Guide.xlabel("Re(θ) [deg]"), Guide.ylabel("Im(θ) [deg]"))

    draw(SVG("naa_atrefht.svg", 13inch, 8inch), p)
end

function prepareplotdata(edges, quadrants, phasediffs, geom2fcn)
    edges = collect(edges)
    numedges = length(edges)
    @assert numedges == length(phasediffs)

    edgecolors = ones(Int, numedges)*9999
    for ii in eachindex(edges)
        if phasediffs[ii] == 2  # suspected edges
            edgecolors[ii] = 5
        elseif phasediffs[ii] == 0
            edgecolors[ii] = quadrants[getindex(geta(edges[ii]))]
        end
    end

    edgecolors = repeat(edgecolors, inner=3)  # x, y is xₐ, xᵦ, NaN, repeat

    x, y = getplotxy(edges)
    z = geom2fcn.(x, y)
    df = DataFrame(x = real(z), y = imag(z), c = edgecolors)

    sort!(df, :c)

    return df
end
