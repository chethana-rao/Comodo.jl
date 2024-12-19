using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Statistics
using LinearAlgebra

# 
Z_cut_offset_top = 200.0
Z_cut_level_skin = 160.0
Z_cut_level_tibia = Z_cut_level_skin+50
Z_cut_level_fibula = Z_cut_level_skin+10

Z_thickness_distal = 20 

fileName_set = ("/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/post/FMA7163_right_leg_isolated_remesh_25k.stl",
"/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/post/FMA24474_10k.stl",
"/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/post/FMA24486_10k.stl",
"/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/post/FMA24477_10k.stl",
"/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/post/FMA24480_10k.stl")
# "/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/BodyParts3D_data/stl/FMA24474.stl",
# "/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/BodyParts3D_data/stl/FMA24486.stl",
# "/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/BodyParts3D_data/stl/FMA24477.stl",
# "/home/simulimb/Che_bodyparts3d/BodyParts3D/assets/BodyParts3D_data/stl/FMA24480.stl")


nameSet = ("skin","femur","patella","tibia","fibula")
indSkin = 1
indFemur = 2
indPatella = 3
indTibia = 4
indFibula = 5

nameSetCut_top = ("skin","femur")
# nameSetCut_bottom = ("skin","femur","tibia","fibula")

n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

# FMAID 	English name
# 7163  	Skin
# 16586 	Right hip bone
# 24474 	Right femur
# 24486     Right patella
# 24477 	Right tibia
# 24480 	Right fibula

FF = []
VV = []
for i in 1:1:length(fileName_set)
    M = load(fileName_set[i])
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    V = [Point{3,Float64}(v) for v in V]
    F,V,_,_ = mergevertices(F,V)
    push!(FF,F)
    push!(VV,V)
end

FF_ori = deepcopy(FF)
VV_ori = deepcopy(VV)

# Cutting surfaces 
v_mid_patella = mean(VV[indPatella])

p = [v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]+Z_cut_offset_top] # Point on cutting plane

for modelName in nameSetCut_top
    indexNow = findall(nameSet.==modelName)[1]
    Fn,Vn,_ = trisurfslice(FF[indexNow],VV[indexNow],n,p; output_type=:below, snapTolerance = snapTolerance) 
    FF[indexNow] = Fn
    VV[indexNow] = Vn
end

p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_level_skin) # Point on cutting plane
Fn,Vn,_ = trisurfslice(FF[indSkin],VV[indSkin],n,p; output_type=:above, snapTolerance = snapTolerance) 
FF[indSkin] = Fn
VV[indSkin] = Vn

p = [v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_level_tibia] # Point on cutting plane
Fn,Vn,_ = trisurfslice(FF[indTibia],VV[indTibia],n,p; output_type=:above, snapTolerance = snapTolerance) 
FF[indTibia] = Fn
VV[indTibia] = Vn

p = [v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_level_fibula] # Point on cutting plane
Fn,Vn,_ = trisurfslice(FF[indFibula],VV[indFibula],n,p; output_type=:above, snapTolerance = snapTolerance) 
FF[indFibula] = Fn
VV[indFibula] = Vn


V_tibia = VV[indTibia]
Eb_tibia = boundaryedges(FF[indTibia])
indCurve_tibia = edges2curve(Eb_tibia)
indCurve_tibia = indCurve_tibia[1:end-1]
v_mid_low_tibia = mean(V_tibia[indCurve_tibia])

p_distal = Point{3,Float64}(v_mid_low_tibia[1],v_mid_low_tibia[2],v_mid_low_tibia[3]-Z_thickness_distal)

zlevel = v_mid_patella[3]

Fs = FF[1]
Vs = VV[1]

Eb = boundaryedges(Fs)
Eb_low = Vector{eltype(Eb)}()
for e in Eb  
       println(Vs[e[1]][3])    
    if Vs[e[1]][3] < zlevel        
      push!(Eb_low,e)  
    end
end
indCurve = edges2curve(Eb_low)
indCurve = indCurve[1:end-1] 

NV = vertexnormal(Fs,Vs; weighting=:area) # May be inefficient to do this to all vertices (/faces)
NV_indCurve = NV[indCurve] # Vertex normals for the curve
NE_indCurve = Vector{Point{3,Float64}}(undef,length(indCurve))
NB_base = Vector{Point{3,Float64}}(undef,length(indCurve))
m = length(indCurve)
for i in eachindex(indCurve)
    ii = indCurve[i]
    jj = indCurve[mod1(i+1,m)]
    NE_indCurve[i] = normalizevector(Vs[jj] - Vs[ii])

    NB_base[i] = cross(NE_indCurve[i],NV_indCurve[i])
end

numPointBezier = 250
v1 = 50
v2 = 25
VV_B = []
Vb = Vector{Point{3,Float64}}()
for i in eachindex(indCurve)
    pb = [Vs[indCurve[i]], Vs[indCurve[i]] + v1*NB_base[i], p_distal + v2*NV_indCurve[i],p_distal]
    println(typeof(pb))
    
    vb = nbezier(pb,numPointBezier)
    push!(VV_B,vb)
    append!(Vb,vb)
end

Fb,Vb = loftpoints2surf(Vb,length(indCurve);close_loop=false,face_type=:quad)
# Vb = simplexcenter(Eb,V)

## Visualization

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Imported mesh")
for i in 1:1:length(fileName_set)
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0.25)
end

# hp2 = wireframe!(ax1,GeometryBasics.Mesh(Vs,Eb_low), linewidth=5,color=:red)


hp2 = lines!(Vs[indCurve],color=:black,linewidth=6)
# hp1 = poly!(ax1,GeometryBasics.Mesh(Vn,Fn), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0.25)

hp3 = dirplot(ax1,Vs[indCurve],NV_indCurve.*mean(edgelengths(Fs,Vs)); color=:red,linewidth=3,scaleval=1.0,style=:from)
hp4 = dirplot(ax1,Vs[indCurve],NE_indCurve.*mean(edgelengths(Fs,Vs)); color=:green,linewidth=3,scaleval=1.0,style=:from)
hp5 = dirplot(ax1,Vs[indCurve],NB_base.*mean(edgelengths(Fs,Vs))*3; color=:blue,linewidth=3,scaleval=1.0,style=:from)

hp6 = lines!(V_tibia[indCurve_tibia],color=:yellow,linewidth=6)
hp7 = mesh!(ax1,GeometryBasics.Mesh(Vb,Fb), color=:blue, shading = FastShading, transparency=false)

# for i in eachindex(indCurve)
#     lines!(VV_B[i],color=:black,linewidth=1)
# end

fig