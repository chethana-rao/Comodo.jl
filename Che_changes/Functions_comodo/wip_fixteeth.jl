using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Comodo.Rotations


testCase = 1
n = 3
r = 5.0

if testCase == 1 
    F,V = geosphere(n,r)
    VF = simplexcenter(F,V)
    B = [vf[1]<r/2.0 for vf in VF] 
elseif testCase == 2
    F,V = quadsphere(n,r)
    Q = RotXYZ(0.0,0.25*π,0.25*π)
    # V = [GeometryBasics.Point{3, Float64}(Q*v) for v ∈ V] 


    VF = simplexcenter(F,V)
    B = [vf[1]<r/2.0 for vf in VF]        
end 


function fixteeth(F,B; method=:add)

    if method == :add
        # Add interior teet 
        indNodesOut = unique(reduce(vcat,F[B]))
        B_fixed = Vector{Bool}(undef,length(B))
        for i in eachindex(B)
            B_fixed[i] = all([in(i,indNodesOut) for i in F[i]])
        end
    elseif method == :remove
        indNodesOut = unique(reduce(vcat,F[.!B]))
        B_fixed = Vector{Bool}(undef,length(B))
        for i in eachindex(B)
            B_fixed[i] = !all([in(i,indNodesOut) for i in F[i]])
        end
    # else
    # error
    end
    return B_fixed
end

B_fix1 = fixteeth(F,B; method=:add)
B_fix2 = fixteeth(F,B; method=:remove)

################ 

## Visualization
cmap = cgrad(:viridis, 2, categorical = true)
strokewidth = 1
linewidth = 3

Fp,Vp = separate_vertices(F,V)
Cp = simplex2vertexdata(Fp,B)
Cp2 = simplex2vertexdata(Fp,B_fix1)
Cp3 = simplex2vertexdata(Fp,B_fix2)

fig = Figure(size=(1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
hp2 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
hp3 = poly!(ax3,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp3, shading = FastShading, colormap=cmap)


fig