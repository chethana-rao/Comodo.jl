using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO


ind_ins= Cb .== 1 
indsh = findall(ind_ins)
Finner = Fbs[indsh]
Finner,Vinner = separate_vertices(Finner,Vbs)            
#try isolating Finner and Vinner 
Ms = GeometryBasics.Mesh(Vinner,Finner)
#visualizatio
f = Figure(size=(100,100))
ax1 = Axis3(f[1, 1],aspect = :data,title="Cylinder")
hp1 = mesh!(ax1,Ms, shading = FastShading, transparency=true, overdraw=false)
display(GLMakie.Screen(),f)


# open("Voutter.csv", "w") do io
#     for v in Vinner
#         println(io, "$(v.data[1]),$(v.data[2]),$(v.data[3])")
#     end
# end

# open("Foutter.csv", "w") do io
#     for f in Finner
#         println(io, "$(f.data[1]),$(f.data[2]),$(f.data[3])")
#     end
# end
