using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Comodo.Rotations


n = 3
r = 5.0

F,V = geosphere(n,r)    
VC = simplexcenter(F,V)
# B = [vf[2]<r/5.0 for vf in VF]
searchTol = r/5.0 # Tolerance for cropping hemisphere (safe, somewhat arbitrary if much smaller then mesh edge lengths)
F = [F[i] for i in findall(map(v-> v[3]>-0.5,VC))]
# F = [F[i] for i in findall(map(f-> mean([V[j][1] for j in f])<=searchTol,F))] # Remove faces below equator
F,V,_ = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed
Eb = boundaryedges(meshedges(F)) # or equivalently Eb = boundaryedges(meshedges(F))
ind = edges2curve(Eb)

B = Vector{Bool}(undef,length(F))
    for i in eachindex(F)
        B[i] =  in(i,ind)
    end





indNodesOut = unique(reduce(vcat,F[B]))
B_fixed = Vector{Bool}(undef,length(B))
    for i in eachindex(F)
        B_fixed[i] = all([in(i,indNodesOut) for i in F[i]])
    end
ind_new = []
for i in eachindex(B_fixed)
    if B_fixed[i] == 1
        push!(ind_new,i)
    end
end
    # Fp = F[B_fixed]

# Fp,Vp = separate_vertices(F,V)
# Cp = simplex2vertexdata(Fp,B)
# Cp2 = simplex2vertexdata(Fp,B_fixed)
# Cp3 = simplex2vertexdata(Fp,B_fix2)
F_new = F[B_fixed]
Eb2 = boundaryedges(F_new) # or equivalently Eb = boundaryedges(meshedges(F))
# ind_new2 = edges2curve(Eb2)

#Visualization
cmap = cgrad(:viridis, 2, categorical = true)
strokewidth = 1
linewidth = 3


fig1 = Figure(size=(1200,800))
ax1 = Axis3(fig1[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)
# ax2 = Axis3(fig1[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
hp3 = scatter!(ax1, V[ind],markersize=15,color=:orange)
ax2 = Axis3(fig1[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
hp2 = poly!(ax2,GeometryBasics.Mesh(V,F_new), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
hp2 = wireframe!(ax2,GeometryBasics.Mesh(V,Eb2), linewidth=5,color=:red)
hp3 = scatter!(ax2, V[ind_new],markersize=15,color=:orange)



# hp3 = poly!(ax2,GeometryBasics.Mesh(Vp,F_ind), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
# hp4 = wireframe!(ax2,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)

fig1