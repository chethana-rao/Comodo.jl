n = 3
r = 5.0

F,V = geosphere(n,r)    
VF = simplexcenter(F,V)
# p = r/1.5
B = [vf[3]<r/1.5 for vf in VF]
F_cut = F[B]
F_cut,V_cut,_ = remove_unused_vertices(F_cut,V)
E_cut = boundaryedges(F_cut)

B_fix1 = fixteeth(F,B; method=:remove)
F_true1 = F[B_fix1]
F_true1,V_true1,_ = remove_unused_vertices(F_true1,V)
E1 = boundaryedges(F_true1)


F2 = F[B_fix1]
F2,V2,_ = remove_unused_vertices(F2,V)
E2 = boundaryedges(F2)
ind1 = edges2curve(E2)
B2 = Vector{Bool}(undef,length(F2))
Vedge = V2[ind1]
# Vedge = [v[1] .- 0.22, v[2] .- 0.22, v[3] for v in Vedge]
# B2 = falses(length(F2))  
for i in eachindex(F2)
    f = F2[i]
    v = V2[f]
   B2[i] = any(in(Vedge),v)
end  
B2 = .!B2
B_fix2 = fixteeth(F2,B2; method=:remove)
F_true2 = F2[B_fix2]
F_true2,V_true2,_ = remove_unused_vertices(F_true2,V2)
E2 = boundaryedges(F_true2)

Fp,Vp = separate_vertices(F,V)
Cp = simplex2vertexdata(Fp,B)
cmap = cgrad(:viridis, 2, categorical = true)
strokewidth = 1
linewidth = 3

Cp2 = simplex2vertexdata(Fp,B_fix1)
# Cp3 = simplex2vertexdata(Fp,B_fix2)



fig = Figure(size=(1200,800))

# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
# hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)
# hp2 = wireframe!(ax1,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)

ax2 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
hp3 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)
hp4 = wireframe!(ax2,GeometryBasics.Mesh(V_cut,E_cut), linewidth=5,color=:red)
hp5 = wireframe!(ax2,GeometryBasics.Mesh(V_true1,E1), linewidth=5,color=:cyan)
# hp5 = scatter!(ax2,V2[ind1], markersize=15,color=:blue)
hp8 = wireframe!(ax2,GeometryBasics.Mesh(V_true2,E2), linewidth=5,color=:blue)
# Legend(fig[2, 1][1,2],[hp1,hp2],["Initial","Refined"])
# Legend(fig[1, 3][4,2], [hp2, hp5, hp8], ["Before","After: Added", "After: Removed"])
display(GLMakie.Screen(),fig)
