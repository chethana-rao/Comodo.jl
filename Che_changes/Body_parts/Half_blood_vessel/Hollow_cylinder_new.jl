using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO


#parameter
pointSpacing = 2.0
r = 10.0 # Radius
n = ceil(Int,(2.0*π*r)/pointSpacing)
height = 5.0*r
r2 = r*0.75
num_steps = ceil(Int,height/pointSpacing)
if iseven(num_steps)
    num_steps += 1
end
close_loop = true
#Outer cylinder_bottom
V1 = circlepoints(r,n;dir=:acw)
r2 = r*0.75
V2 = circlepoints(r2,n;dir=:acw)
VT = (V1,V2) # Curves
R = ([1,2],) # hollow
# R =([1,2],[2]) #two regions
P = (pointSpacing,pointSpacing)  # Point spacings
Fb,Vb,Cb = regiontrimesh(VT,R,P)
invert_faces!(Fb)
#Outer cylinder top
V1_t = circlepoints(r,n;dir=:acw)
V1_t = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V1_t]
V2_t = circlepoints(r2,n;dir=:acw)
V2_t = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V2_t] 
VT = (V1_t,V2_t) # Curves
R = ([1,2],) # hollow
P = (pointSpacing,pointSpacing)  # Point spacings
Ft,Vt,Ct = regiontrimesh(VT,R,P)
#loft linear 
Fl_1,Vl_1 = loftlinear(V1,V1_t;num_steps=num_steps,close_loop=close_loop,face_type =:tri)
Fl_2,Vl_2 = loftlinear(V2,V2_t;num_steps=num_steps,close_loop=close_loop,face_type =:tri)
Fhc,Vhc,Chc = joingeom(Fb,Vb,Ft,Vt,Fl_1,Vl_1,Fl_2,Vl_2)
Fhc,Vhc,Chc = joingeom(Fb,Vb,Ft,Vt,Fl_1,Vl_1)

Fhc,Vhc = mergevertices(Fhc,Vhc)
D1 = edgelengths(Fhc,Vhc)
#inner cylinder 
Vo = circlepoints(r,n;dir=:acw)
V1 = circlepoints(r2,n;dir=:acw)
Vr2 = circlepoints(r2,n;dir=:acw)
V2 = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V2]
VT = (Vo,V1,)
R = ([2],)
P = (pointSpacing)
Fb1,Vb1,_ = regiontrimesh(VT,R,P)
invert_faces!(Fb1)
Vo = circlepoints(r,n;dir=:acw)
Vo = [GeometryBasics.Point{3, Float64}(v[1],v[2],height*0.9) for v ∈ Vo]
VT = (Vo,V2,)
R = ([2],)
P = (pointSpacing)
Ft,Vt,_ = regiontrimesh(VT,R,P)
# Fnew,Vnew,Cnew = joingeom()
Fic,Vic,Cic = joingeom(Fb1,Vb1,Fl_2,Vl_2,Ft,Vt)
Fic,Vic = mergevertices(Fic,Vic)
D2 = edgelengths(Fic,Vic)

#visualization
Fw,Vw,Cw = joingeom(Fhc,Vhc,Fic,Vic)
Fw,Vw = mergevertices(Fw,Vw)

# Visualisation hollow cylinder

fig = Figure(size=(1200,1000))
ax1 = Axis3(fig[1, 1],aspect = :data,title="Cylinder")
ax2 = Axis3(fig[1, 2],aspect = :data,title="Merged")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vhc,Fhc), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp2 = poly!(ax1,GeometryBasics.Mesh(Vic,Fic), strokewidth=1,color=:green, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp3 = poly!(ax2,GeometryBasics.Mesh(Vw,Fw), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
normalplot(ax2,Fw,Vw)
display(GLMakie.Screen(),fig)

#
#visualization
# Fic = [f.+length(Vhc) for f in invert_faces(Fic)]    
# Fw = [Fhc;Fic]
#  Vw = [Vhc;Vic]

# Cw = ones(length(Fw))
# Cw[length(Fhc)+1:end] .+= 1
stringOpt = "paAqY"
v_region1 = Point{3,Float64}(r2+((r-r2)/2), 0.0, height/2)
v_region2 = Point{3,Float64}(0.0, 0.0, height/2)
vol1 = mean(D1)^3 / (6.0*sqrt(2.0)) #volume of ideal tetrahedron
vol2 = mean(D2)^3 / (6.0*sqrt(2.0))
region_vol = [vol1,vol2]
E,V,CE,Fb,Cb= tetgenmesh(Fw,Vw; facetmarkerlist=Cw, V_regions=[v_region1,v_region2],region_vol=region_vol, stringOpt)
# E1,V2,CE2,Fb2,Cb2= tetgenmesh(Fic,Vic; facetmarkerlist=Cic, V_regions=[v_region2],region_vol=region_vol[2], stringOpt)

# #Visualization
# ## Visualization tetgen 
cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # Triangular faces
CE_F = repeat(CE,inner=4) #why is it repeating

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)

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
ZE = [v[1] for v in VE]
Z = [v[1] for v in V]
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


#  fileName = comododir()*"/assets/img/TetGen_cylinder_example_01.mp4"
#  slider2anim(fig3,hSlider,fileName; backforth=true, duration=6)
