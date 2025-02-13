using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Comodo.Rotations


n = 3
r = 5.0

F,V = geosphere(n,r)    
VF = simplexcenter(F,V)
B = [vf[3]<r/1.5 for vf in VF]
F_cut = F[B]
F_cut,V_cut,_ = remove_unused_vertices(F_cut,V)
E_cut = boundaryedges(F_cut)
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


F_true1 = F[B_fix1]
F_true1,V_true1,_ = remove_unused_vertices(F_true1,V)
E1 = boundaryedges(F_true1)

F_true2 = F[B_fix2]
F_true2,V_true2,_ = remove_unused_vertices(F_true2,V)
E2 = boundaryedges(F_true2)
Fp,Vp = separate_vertices(F,V)
Cp = simplex2vertexdata(Fp,B)
#Visualization
cmap = cgrad(:viridis, 2, categorical = true)
strokewidth = 1
linewidth = 3

Cp2 = simplex2vertexdata(Fp,B_fix1)
Cp3 = simplex2vertexdata(Fp,B_fix2)

fig = Figure(size=(1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)
hp8 = wireframe!(ax1,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
hp2 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)
hp3 = wireframe!(ax2,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)
hp4 = wireframe!(ax2,GeometryBasics.Mesh(V_true1,E1), linewidth=5,color=:cyan)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
hp5 = poly!(ax3,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp3, shading = FastShading, colormap=cmap)
hp6 = wireframe!(ax3,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)
hp7 = wireframe!(ax3,GeometryBasics.Mesh(V_true2,E2), linewidth=5,color=:cyan)

fig

# ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
# hp4 = poly!(ax4,GeometryBasics.Mesh(V_cut,F_cut), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
# hp5 = wireframe!(ax4,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)

# ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
# hp6 = poly!(ax5,GeometryBasics.Mesh(V_true1,F_true1), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
# hp7 = wireframe!(ax5,GeometryBasics.Mesh(V_true1,E1), linewidth=5,color=:red)

# ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
# hp8 = poly!(ax6,GeometryBasics.Mesh(V_true2,F_true2), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
# hp9 = wireframe!(ax6,GeometryBasics.Mesh(V_true2,E2), linewidth=5,color=:red)


fig





#Visualization


# fig1 = Figure(size=(1200,800))
# ax1 = Axis3(fig1[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
# hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color =Cp, shading = FastShading, colormap=cmap)
# # hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)
# # ax2 = Axis3(fig1[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
# # hp3 = scatter!(ax1, V[ind],markersize=15,color=:orange)
#  ax2 = Axis3(fig1[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
#  hp2 = poly!(ax2,GeometryBasics.Mesh(V_true,F_true), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
# # hp2 = wireframe!(ax2,GeometryBasics.Mesh(V,Eb2), linewidth=5,color=:red)
# # hp3 = scatter!(ax2, Vp[ind_false],markersize=15,color=:orange)



# # hp3 = poly!(ax2,GeometryBasics.Mesh(Vp,F_ind), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
#  hp4 = wireframe!(ax2,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)

# fig1