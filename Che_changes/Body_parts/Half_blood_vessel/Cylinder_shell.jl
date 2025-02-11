using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO

"""
    cylinder_shell(n::T,r::T, h::T; face_type=:tri,close_loop, closed=false) where T <: Real

Lofts surface between curves

# Description 

The `cylinder_shell` function utilises a `loftlinear`operation to create a cylinder surface. If `closed==true`, 
then it is assumed the top and botom surafces of the cylinder are closed, and `regiontrimesh` is used for closing 
the cylinder. `closed == false` implies a simple extrude of the circle. The user can request different face types 
for the cylinder surface. The default is `face_type=:tri` which will form isoceles triangles (or equilateral 
triangles if the spacing is even) for a planar curve. The other `face_type` options supported are `:quad` (quadrilateral),
and `:tri_slash`. For the latter, triangles are formed by slashing the quads.  

# Arguments:
- `n`:: T pointspacing
- `r:: T` Radius 
- `h :: T`Height 
"""


function cylinder_shell(n::T,r::T, h::T; face_type=:tri,close_loop, closed=false) where T <: Real
    if n<0
        throw(ArgumentError("n should be >= 0"))
    end
    n1 = ceil(Int,(2.0*π*r)/n)
    V1 = circlepoints(r,n1;dir=:acw)
    V2 = circlepoints(r,n1;dir=:acw)
    V2 = [GeometryBasics.Point{3, Float64}(v[1],v[2],h) for v ∈ V1]
    num_steps = ceil(Int,h/n) + 1
    if iseven(num_steps)
        num_steps += 1
    end
    # close_loop = true
    if closed == true      
    Fl,Vl = loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type =face_type)
    VT = (V1,)
    R = ([1],)
    P = (n)
    Fb,Vb,_ = regiontrimesh(VT,R,P)
    invert_faces!(Fb)
    VT= (V2,)
    R = ([1],)
    P = (n)
    Ft,Vt,_ = regiontrimesh(VT,R,P)  
    Fc,Vc,Cc = joingeom(Fb,Vb,Fl,Vl,Ft,Vt)
    Fc,Vc = mergevertices(Fc,Vc)
    elseif closed == false 
        Fc,Vc = loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type =face_type)
        Cc = ones(length(Vc))
    else throw(ArgumentError("Closed should be either true of false "))
    end 
    return Fc,Vc,Cc
end 
#demo
n = 2.0
r = 5.0 # Radius   
h = 5*r
face_types = [:tri, :tri_even, :backslash, :forwardslash, :quad, :quad2tri]
close_loop = true
testcase = 2
if testcase == 1 
    closed = true
elseif testcase ==2
     closed = false
end 
     #Visualisation
markersize = 10
linewidth = 2 

Fig = Figure(size=(1200,1200))

nRows = 2
for (q,face_type) in enumerate(face_types)
    i = mod1(q,nRows)
    j = ceil.(Int,q/nRows)                
    F,V,_ =cylinder_shell(n,r, h; face_type= face_type,close_loop = close_loop, closed=closed) 
    ax1 = Axis3(Fig[i, j], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Cylinder Surface, Covered:$closed, face_type=:$face_type")    
    hp=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white,shading=FastShading,transparency=false)
    # normalplot(ax1,GeometryBasics.Mesh(V,F); type_flag=:face, color=:black)
    j+=1
end
# Legend(fig[1, 4],[hp1,hp2,hp3],["curve 1", "curve 2", "lofted surface"])

display(GLMakie.Screen(),Fig)
# t = 2
# E, Ve = extrudefaces(F,V; extent=t, direction=:positive, num_steps=n)


