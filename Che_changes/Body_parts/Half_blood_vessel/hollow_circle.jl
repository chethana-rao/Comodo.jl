using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO



pointSpacing = 2.0
r = 10.0 # Radius
n = ceil(Int,(2.0*π*r)/pointSpacing)
V1 = circlepoints(r,n;dir=:acw)
r2 = r*0.75
V2 = circlepoints(r2,n;dir=:acw)
VT = (V1,V2) # Curves
R = ([1,2],) # hollow
# R =([1,2],[2]) #two regions
P = (pointSpacing,pointSpacing)  # Point spacings
Fb,Vb,Cb = regiontrimesh(VT,R,P)
invert_faces!(Fb)
#top
height = 5.0*r
V1_t = circlepoints(r,n;dir=:acw)
V1_t = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V1_t]

r2 = r*0.75
V2_t = circlepoints(r2,n;dir=:acw)
V2_t = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V2_t] 
VT = (V1_t,V2_t) # Curves
R = ([1,2],) # hollow
P = (pointSpacing,pointSpacing)  # Point spacings
Ft,Vt,Ct = regiontrimesh(VT,R,P)
#loft linear 
num_steps = ceil(Int,height/pointSpacing)
if iseven(num_steps)
    num_steps += 1
end
close_loop = true
Fl_1,Vl_1 = loftlinear(V1,V1_t;num_steps=num_steps,close_loop=close_loop,face_type =:tri)
Fl_2,Vl_2 = loftlinear(V2,V2_t;num_steps=num_steps,close_loop=close_loop,face_type =:tri)
Fhc,Vhc,Chc = joingeom(Fb,Vb,Ft,Vt,Fl_1,Vl_1,Fl_2,Vl_2)
Fhc,Vhc = mergevertices(Fhc,Vhc)


# Visualisation hollow cylinder

fig = Figure(size=(1200,1000))
ax1 = Axis3(fig[1, 1],aspect = :data,title="Hollow region")
ax2 = Axis3(fig[1, 2],aspect = :data,title="merged")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vt,Ft), strokewidth=1,color=:blue, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp2 = poly!(ax1,GeometryBasics.Mesh(Vb,Fb), strokewidth=1,color=:green, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp3 = poly!(ax1,GeometryBasics.Mesh(Vl_1,Fl_1), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp4 = poly!(ax1,GeometryBasics.Mesh(Vl_2,Fl_2), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp4 = poly!(ax2,GeometryBasics.Mesh(Vhc,Fhc), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
normalplot(ax2,Fhc,Vhc)
display(GLMakie.Screen(),fig)

#tetgen
stringOpt = "paAqY"
v_region = Point{3,Float64}(r2+((r-r2)/2), 0.0, height/2)
vol1 = pointSpacing^3 / (6.0*sqrt(2.0)) #volume of ideal tetrahedron
E,V,CE,Fb1,Cb1= tetgenmesh(Fhc,Vhc; facetmarkerlist=Chc, V_regions=[v_region],region_vol=vol1, stringOpt)
## Visualization tetgen 
cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # Triangular faces
CE_F = repeat(CE,inner=4)

Fbs,Vbs = separate_vertices(Fb1,V)
Cbs_V = simplex2vertexdata(Fbs,Cb1)

Fs,Vs = separate_vertices(F,V)
CE_Vs = simplex2vertexdata(Fs,CE_F)
M = GeometryBasics.Mesh(Vs,Fs)

strokewidth = 1 


fig3 = Figure(size=(800,800))

ax1 = Axis3(fig3[1, 1][1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary surfaces")
hp1 = mesh!(ax1,GeometryBasics.Mesh(Vbs,Fbs), color=Cbs_V, shading = FastShading, transparency=true, overdraw=false,colorrange = (1,3),colormap=cmap)
# scatter!(ax1,V_regions,color=:yellow,markersize=25)
# scatter!(ax1,V_holes,color=:red,markersize=25)

ax2 = Axis3(fig3[1, 1][1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Cut mesh")
hp2 = poly!(ax2,M, color=CE_Vs, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=strokewidth, overdraw=false,colorrange = (1,2),colormap=cmap)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F,V)))

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider = Slider(fig3[2, 1], range = stepRange, startvalue = mean(stepRange),linewidth=30)

on(hSlider.value) do z 

    B = ZE .<= z
    indShow = findall(B)
    if isempty(indShow)
        hp2.visible=false        
    else        
        hp2.visible=true
        Fs = element2faces(E[indShow])
        Cs = repeat(CE[indShow],inner=4)
        
        indB = boundaryfaceindices(Fs)        
        Fs = Fs[indB]
        Cs = Cs[indB]
        Fs,Vs = separate_vertices(Fs,V)
        CE_Vs = simplex2vertexdata(Fs,Cs)
        Ms = GeometryBasics.Mesh(Vs,Fs)
        hp2[1] = Ms
        hp2.color = CE_Vs
    end

end
# hSlider.selected_index[]+=1
slidercontrol(hSlider,ax2)

display(GLMakie.Screen(),fig3)

