module RolyVis

import Luxor
using Plots: @colorant_str, palette, color, Color
using StaticArrays
using Statistics: mean

using Roly: AssemblySystem, Polyform, species, nsites
using Roly
using Graphs
using LinearAlgebra

export DEFAULT_COLORS, _draw_polyform, draw_polyform_grid, draw_polyforms

const DEFAULT_COLORS = [
    [colorant"#B2D0EA", colorant"#54A3E4", colorant"#1A78C6", colorant"#00549A"],
    [colorant"#F39F9D", colorant"#ED7E7C", colorant"#E75451", colorant"#BE2E2C"],
    [colorant"#FEE3B5", colorant"#FFC563", colorant"#FFB12F", colorant"#FFA000"],
    [colorant"#DBB9E4", colorant"#D99DE8", colorant"#B571C4", colorant"#8E4E9C"],
    [colorant"#AAF5FF", colorant"#86E4F1", colorant"#43CFE2", colorant"#20BCD2"],
    [colorant"#BCF3E1", colorant"#3DD4A4", colorant"#19B684", colorant"#008D60"],
    [colorant"#FFD4B6", colorant"#F7AF7D", colorant"#FF8835", colorant"#ED6304"],
    [colorant"#C5C9FE", colorant"#A9AFFF", colorant"#7C85FF", colorant"#4854FD"],
    [colorant"#F4FFC4", colorant"#E0F290", colorant"#B9D63E", colorant"#859F13"],
    [colorant"#FFB2C5", colorant"#FF8EA9", colorant"#FF5C83", colorant"#E32654"],
    [colorant"#95FABA", colorant"#4DD980", colorant"#21AB53", colorant"#008832"]
]

function scaled_colors(base_color, n)
    return base_color .* [1 - i / n for i in 0:(n - 1)]
end

function draw_ngon(n, radius; colors="grey90", text=nothing, textcolor=colorant"black", holescale=5//8, cornerscale=30, linescale=30)
    if colors isa Color
        colors = fill(color, n)
    end

    a_out = radius
    a_in = radius * holescale
    corner_r = radius / cornerscale * n

    angles = [i * 2π / n for i in 0:(n - 1)]
    orientation = 2π / n * (floor(n // 2) - 1) - π / 2
    if iseven(n)
        orientation -= π / n
    end

    interior_angle = (n - 2) / n * π

    Luxor.setline(radius / linescale)

    Luxor.@layer begin
        for (c, θ) in zip(colors, angles)
            Luxor.@layer begin
                Luxor.rotate(θ)
                ngon_outer = Luxor.ngonside(Luxor.O, a_out, n, orientation; vertices=true)
                clip_points = 1.2 * [ngon_outer[2], Luxor.Point(0, 0), ngon_outer[1]]
                Luxor.poly(clip_points, :clip)

                Luxor.sethue(c)
                Luxor.polysmooth(ngon_outer, corner_r, :fill)
            end
        end

        # Draw the outline
        Luxor.sethue("black")
        ngon_outer = Luxor.ngonside(Luxor.O, a_out, n, orientation; vertices=true)
        Luxor.polysmooth(ngon_outer, corner_r, :stroke)

        # Draw the hole in the center
        ngon_inner = Luxor.ngonside(Luxor.O, a_in, n, orientation; vertices=true)

        Luxor.sethue("white")
        Luxor.polysmooth(ngon_inner, corner_r, :fill)

        Luxor.sethue("black")
        Luxor.polysmooth(ngon_inner, corner_r, :stroke)

        # Draw lines delimiting the binding sides
        offset = 1 / sin(interior_angle / 2) - 1
        multiplyer = n == 3 ? 2 : 1 # don't ask me why this is needed
        for (ngi, ngo) in zip(ngon_inner, ngon_outer)
            Luxor.line(ngi * (1 - multiplyer * corner_r / a_in * offset),
                       ngo * (1 - multiplyer * corner_r / a_out * offset), :stroke)
        end

        Luxor.fontsize(radius/2)
        isnothing(text) || Luxor.text(string(text), Luxor.O + Luxor.Point(0, radius/6), halign=:center, valign=:center)
    end
end

function _draw_polyform(pform, assembly_system; size=200, colors=nothing, label_species=false)
    spes = species(pform)
    np, nb = Roly.size(assembly_system)
    geoms = Roly.geometries(assembly_system)

    if isnothing(colors)
        if np <= length(DEFAULT_COLORS)
            colors = DEFAULT_COLORS
        elseif np < 20
            colors = palette(:tab20)
        else
            colors = palette(:plasma, np)
        end
    end

    d = size
    xs = [d * SVector(1, -1) .* x for x in pform.xs]
    ψs = [ψ.θ for ψ in pform.ψs]

    polygons = [nsites(assembly_system.geometries[k])
                for k in species(pform)]

    # Make sure the center of mass is centered
    x_com = mean(xs)
    xs .-= Ref(x_com)

    # Draw all the triangles at the correct positions
    for (particle, (x, ψ)) in enumerate(zip(xs, ψs))
        n = nsites(geoms[particle])
        base_colors = colors[spes[particle]]
        if base_colors isa Color
            base_colors = scaled_colors(base_colors, n)
        end
        color = [Roly.isinert(spes[particle], i, assembly_system) ? colorant"#E7E7E7" : base_colors[i] for i in 1:n]

        Luxor.@layer begin
            Luxor.translate(x...)
            Luxor.rotate(-ψ * π)

            draw_ngon(polygons[particle], size; colors=color, text=label_species ? spes[particle] : nothing)
        end
    end

    es = edges(pform.anatomy)
    # pedges = [e for e in es if reverse(e) in es && e.src <= e.dst]
    
    # geoms = Roly.geometries(assembly_system)
    # i = 1
    # for edge in pedges
    #     pos_i = Roly.get_sitepos(pform, assembly_system, edge.src)
    #     pos_j = Roly.get_sitepos(pform, assembly_system, edge.dst)
    #     xi = Luxor.Point((d * SVector(1, -1) .* pos_i - x_com)...)
    #     xj = Luxor.Point((d * SVector(1, -1) .* pos_j - x_com)...)

    #     if norm(pos_i - pos_j) > 1e-6
    #         Luxor.@layer begin
    #             Luxor.sethue(species_colors[i])
    #             Luxor.setopacity(0.7)

    #             # Luxor.translate(xi...)
    #             Luxor.circle(xi, 8, action = :fill)
    #             # Luxor.translate(xj...)
    #             Luxor.circle(xj, 8, action = :fill)
    #         end
    #         i += 1
    #     end
    # end

    return
end

function draw_polyform_grid(structures, assembly_system; size=50, box_size=5, ncols=4, colors=nothing,
                         structure_names=nothing, label_species=false)
    nstructures = length(structures)
    nrows = nstructures ÷ ncols
    if nstructures % ncols != 0
        nrows += 1
    end

    structures = [structures[i] for i in sort(collect(keys(structures)))]

    tiles = Luxor.Tiler(box_size * size * ncols, box_size * size * nrows, nrows, ncols)

    if isnothing(structure_names)
        i = 1
        for (structure, (pos, n)) in zip(structures, tiles)
            Luxor.@layer begin
                Luxor.translate(pos)
                _draw_polyform(structure, assembly_system; size, colors, label_species)
            end
            i += 1
        end
    else
        i = 1
        for (structure, (pos, n), name) in zip(structures, tiles, structure_names)
            Luxor.@layer begin
                Luxor.translate(pos)
                _draw_polyform(structure, assembly_system; size, colors, label_species)
                Luxor.fontsize(size / 2)
                Luxor.text(string(name), Luxor.Point(-2 * size / 2, -3 * size / 2))
            end
            i += 1
        end
    end
end

function draw_polyform(polyform::Polyform, assembly_system::AssemblySystem, filename="temp.png"; scale=300, colors=nothing, size=(1000, 1000), label_species=false)
    extension = splitext(filename)[end]
    if extension == ".png"
        return Luxor.@png(_draw_polyform(polyform, assembly_system; size=scale, colors, label_species), size[1], size[2], filename)
    elseif extension == ".pdf"
        return Luxor.@pdf(_draw_polyform(polyform, assembly_system; size=scale, colors, label_species), size[1], size[2], filename)
    else
        error("Can only render png or pdf.")
    end
end

function draw_polyforms(polyforms::AbstractVector{<:Polyform}, assembly_system::AssemblySystem, filename="temp.png"; scale=50, ncols=10, box_size=5, colors=nothing, size=(1000, 1000), label_species=false)
    extension = splitext(filename)[end]
    if extension == ".png"
        return Luxor.@png(draw_polyform_grid(polyforms, assembly_system; size=scale, box_size=box_size, ncols=ncols, colors, label_species), size[1], size[2], filename)
    elseif extension == ".pdf"
        return Luxor.@pdf(draw_polyform_grid(polyforms, assembly_system; size=scale, box_size=box_size, ncols=ncols, colors, label_species), size[1], size[2], filename)
    else
        error("Can only render png or pdf.")
    end
end
end # module RolyVis
