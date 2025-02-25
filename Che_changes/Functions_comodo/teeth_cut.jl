using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Comodo.Rotations


"""
    fixteeth(F,B; method=:add)

# Description 

    Removes or addes hanging nodes (unique nodes that don't share neighbouring nodes at the boundary edges). This function is used to adjust the nodes 
when a mesh is sliced or cut. Depending on the users application hanging nodes are either removed or added to make the cut 
more uniform. The input variables are Faces and Boolean vector B, which indicates the faces that are to be removed and those 
that need to be kept. The output is B_fixed, which again is a Boolean vector that fixes the output nodes. 
# Arguments:
- `F`:: Faces
- `B:: Boolean Vector
- `method :: either add or remove 
"""
function fixteeth(F::Array{NgonFace{N,TF}},B::Vector{Bool}; method=:add) where N where TF<:Integer

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
    else(throw(ArgumentError("Method should be :either add or :remove ")))
    
    end
    return B_fixed
end

testCase = 2

if testCase == 1 
    n = 3
    r = 5.0
    
    F,V = geosphere(n,r)    
    VF = simplexcenter(F,V)
    B = [vf[3]<r/1.5 for vf in VF]
    F_cut = F[B]
    F_cut,V_cut,_ = remove_unused_vertices(F_cut,V)
    E_cut = boundaryedges(F_cut)
    
    B_fix1 = fixteeth(F,B; method=:add)
    B_fix2 = fixteeth(F,B; method=:remove)


    F_true1 = F[B_fix1]
    F_true1,V_true1,_ = remove_unused_vertices(F_true1,V)
    E1 = boundaryedges(F_true1)

    F_true2 = F[B_fix2]
    F_true2,V_true2,_ = remove_unused_vertices(F_true2,V)
    E2 = boundaryedges(F_true2)
    Fp,Vp = separate_vertices(F,V)
    Cp = simplex2vertexdata(Fp,B)
    cmap = cgrad(:viridis, 2, categorical = true)
    strokewidth = 1
    linewidth = 3

    Cp2 = simplex2vertexdata(Fp,B_fix1)
    Cp3 = simplex2vertexdata(Fp,B_fix2)

    fig = Figure(size=(1200,800))

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
    hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)
    hp2 = wireframe!(ax1,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)

    ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
    hp3 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)
    hp4 = wireframe!(ax2,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)
    hp5 = wireframe!(ax2,GeometryBasics.Mesh(V_true1,E1), linewidth=5,color=:cyan)

    ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
    hp6 = poly!(ax3,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp3, shading = FastShading, colormap=cmap)
    hp7 = wireframe!(ax3,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)
    hp8 = wireframe!(ax3,GeometryBasics.Mesh(V_true2,E2), linewidth=5,color=:blue)
    Legend(fig[1, 4], [hp2, hp5, hp8], ["Before","After: added", "After: removed"])
    display(GLMakie.Screen(),fig)

elseif testCase ==2 
    n = 3
    r = 5.0
    F1,V1 = geosphere(n,r)
    E,V,CE,Fb,Cb = tetgenmesh(F1,V1) 
    F = element2faces(E) # Triangular faces   
    VF = simplexcenter(F,V)    
    B = [vf[3]<r/1.5 for vf in VF]
    F_cut = F[B]
    F_cut,V_cut,_ = remove_unused_vertices(F_cut,V)
    E_cut = boundaryedges(F_cut)
    B_fix1 = fixteeth(F,B; method=:add)
    B_fix2 = fixteeth(F,B; method=:remove)
    F_true1 = F[B_fix1]
    F_true1,V_true1,_ = remove_unused_vertices(F_true1,V)
    E1 = boundaryedges(F_true1)
    ind1 = edges2curve(E1)
    F_true2 = F[B_fix2]
    F_true2,V_true2,_ = remove_unused_vertices(F_true2,V)
    E2 = boundaryedges(F_true2)
    Fp,Vp = separate_vertices(F,V)
    Cp = simplex2vertexdata(Fp,B)
    #Visualizatio
    cmap = cgrad(:viridis, 2, categorical = true)
    strokewidth = 1
    linewidth = 3

    Cp2 = simplex2vertexdata(Fp,B_fix1)
    Cp3 = simplex2vertexdata(Fp,B_fix2)

    fig = Figure(size=(1200,800))

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
    hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)
    ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
    hp2 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)
    ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
    hp5 = poly!(ax3,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp3, shading = FastShading, colormap=cmap)
    ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
    hp4 = poly!(ax4,GeometryBasics.Mesh(V_cut,F_cut), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
    ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
    hp6 = poly!(ax5,GeometryBasics.Mesh(V_true1,F_true1), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
    ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
    hp8 = poly!(ax6,GeometryBasics.Mesh(V_true2,F_true2), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
    display(GLMakie.Screen(),fig)

end 

