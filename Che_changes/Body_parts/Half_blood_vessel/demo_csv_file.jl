using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO
using Printf


# Finner = Fb[Cb .== 1]
# # Finner,Vinner = separate_vertices(Finner,Vbs)            
# #try isolating Finner and Vinner 
# Ms = GeometryBasics.Mesh(V,Finner)
# #visualizatio
# f = Figure(size=(100,100))
# ax1 = Axis3(f[1, 1],aspect = :data,title="Cylinder")
# hp1 = mesh!(ax1,Ms, shading = FastShading, transparency=true, overdraw=false)
# display(GLMakie.Screen(),f)


open("Nodes_data.csv", "w") do io
    for v in V
         @printf(io,"%.16e, %.16e, %.16e\n",v[1],v[2], v[3])
     #   println(io, "$(v.data[1]),$(v.data[2]),$(v.data[3])")
    end
end

open("Faces_boundary.csv", "w") do io
    for (i,f) in enumerate(Fb)
        @printf(io,"%i, %i, %i, %i\n",f[1],f[2],f[3],Cb[i])
      #  println(io, "$(f.data[1]),$(f.data[2]),$(f.data[3])")
    end
end
open("Elements.csv", "w") do io
    for (i,e) in enumerate(E)
        @printf(io,"%i, %i, %i, %i\n",e[1],e[2],e[3],CE[i])
      #  println(io, "$(f.data[1]),$(f.data[2]),$(f.data[3])")
    end
end
