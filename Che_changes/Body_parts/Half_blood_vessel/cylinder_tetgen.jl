using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO


#create a cylinder using extrude face 
pointSpacing = 2.0
r = 10.0 # Radius

height = 5.0*r
n = ceil(Int,(2.0*π*r)/pointSpacing)
V1 = circlepoints(r,n;dir=:acw)
V2 = circlepoints(r,n;dir=:acw)
V2 = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V2]
num_steps = ceil(Int,height/pointSpacing)
if iseven(num_steps)
    num_steps += 1
end
close_loop = true
Fl,Vl = loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type =:tri)

VT = (V1,)
R = ([1],)
P = (pointSpacing)
Fb,Vb,_ = regiontrimesh(VT,R,P)
invert_faces!(Fb)
VT = (V2,)
R = ([1],)
P = (pointSpacing)
Ft,Vt,_ = regiontrimesh(VT,R,P)

Fc,Vc,Cc = joingeom(Fb,Vb,Fl,Vl,Ft,Vt)
Fc,Vc = mergevertices(Fc,Vc)

#TetGen
stringOpt = "paAqY"
v_region = Point{3,Float64}(0.0, 0.0, height/2)
vol1 = pointSpacing^3 / (6.0*sqrt(2.0)) #volume of ideal tetrahedron
E,V,CE,Fb1,Cb1= tetgenmesh(Fc,Vc; facetmarkerlist=Cc, V_regions=[v_region],region_vol=vol1, stringOpt)


#Visualization
Fp,Vp = separate_vertices(Fc,Vc) # Give each face its own point set 
Cp = simplex2vertexdata(Fp,Cc) # Convert face color data to vertex color data 

fig2 = Figure(size=(1200,1200))
ax1 = Axis3(fig2[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loft linear")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vl,Fl), strokewidth=1,color=:white,shading=FastShading,transparency=false)
hp2 = poly!(ax1,GeometryBasics.Mesh(Vb,Fb), strokewidth=1,color=:green,shading=FastShading,transparency=false)
hp3 = poly!(ax1,GeometryBasics.Mesh(Vt,Ft), strokewidth=1,color=:blue,shading=FastShading,transparency=false)
ax2 = Axis3(fig2[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Merged")
hp1 = poly!(ax2,GeometryBasics.Mesh(Vc,Fc), strokewidth=1,color=:white,shading=FastShading,transparency=false)
# normalplot(ax2,Fc,Vc)
display(GLMakie.Screen(),fig2)
fig = Figure(size=(1200,1000))
ax1 = Axis3(fig[1, 1],aspect = :data,title="Multi-region meshing")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth=1,color=Cp, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
Colorbar(fig[1, 1][1, 2], hp1)

display(GLMakie.Screen(),fig)

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

fig3
