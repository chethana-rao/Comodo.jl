using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

#mesh parameters
nx = 7
num_steps = nx +1
V_quad = gridpoints(range(0.0,1.0,num_steps),range(0.0,1.0,num_steps),0.0)
V_tri = deepcopy(V_quad)
periodicity1 = (false, false)
periodicity2 = (false, false)    
face_type1 = :quad
face_type2 = :tri
F_quad = grid2surf(V_quad,num_steps; face_type=:quad, periodicity=periodicity1)
F_tri= grid2surf(V_tri,num_steps; face_type=:tri, periodicity=periodicity2,tri_dir=1)
#boundary edges
E_quad = boundaryedges(F_quad)
E_tri = boundaryedges(F_tri)
M_quad = GeometryBasics.Mesh(V_quad,F_quad)
M_tri = GeometryBasics.Mesh(V_tri,F_tri)
VG = deepcopy(V_quad)
# Boundary edges 
bottom_elements= E_quad[1:nx]
bottom_nodes = V_quad[elements2indices(bottom_elements)]
right_elements = E_quad[nx+1:2*nx]
right_nodes = V_quad[elements2indices(right_elements)]
top_elements = E_quad[(2*nx+1):3*nx]
top_nodes = V_quad[elements2indices(top_elements)]
left_elements = E_quad[(3*nx+1):4*nx]
left_nodes = V_quad[elements2indices(left_elements)]
corner_indices = [1,nx+1,(nx+1)^2-nx,(nx+1)^2]
corner_nodes = V_quad[corner_indices]
#visualisations
markerSize = 20
strokeWidth = 2
linewidth = 4
fig1 = Figure(size=(1600,1200))
ax1 = Axis3(fig1[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=quad",azimuth=-pi/2,elevation=pi/2)
poly!(ax1,M_quad, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
hp1=wireframe!(ax1,GeometryBasics.Mesh(V_quad,bottom_elements),linewidth=linewidth, transparency=false, color=:red)
hp2=wireframe!(ax1,GeometryBasics.Mesh(V_quad,right_elements),linewidth=linewidth, transparency=false, color=:blue)
hp3=wireframe!(ax1,GeometryBasics.Mesh(V_quad,top_elements),linewidth=linewidth, transparency=false, color=:black)
hp4=wireframe!(ax1,GeometryBasics.Mesh(V_quad,left_elements),linewidth=linewidth, transparency=false, color=:yellow)
# scatter!(ax1,VG,color=:green,markersize=markerSize)
scatter!(ax1,bottom_nodes,color=:black,markersize=markerSize)
scatter!(ax1,right_nodes,color=:green,markersize=markerSize)
scatter!(ax1,top_nodes,color=:yellow,markersize=markerSize)
scatter!(ax1,left_nodes,color=:red,markersize=markerSize)

Legend(fig1[:,2],[hp1,hp2,hp3,hp4],["Bottom","Right","Top","Left"])
ax2 = Axis3(fig1[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=tri",azimuth=-pi/2,elevation=pi/2)
poly!(ax2,M_tri, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
wireframe!(ax2,GeometryBasics.Mesh(V_tri,E_tri),linewidth=linewidth, transparency=false, color=:red)
scatter!(ax2,VG,color=:green,markersize=markerSize)
display(GLMakie.Screen(),fig1)
#Finite Element analysis    
#compute FEA control parameters
V = deepcopy(V_quad)
numNodes = length(V)
NumDof = 2
numSpatialDim = length(V[1])-1
numNodesElem = length(F_quad[1])
numElements = length(F_quad)
numDOF_el = NumDof * numNodesElem
numEq = numNodes * NumDof
#Gaussian Quadrature 
s = [-1 1 1 -1]
t = [-1 -1 1 1]
order = 2
if order == 1
    xi = 0
    eta = 0
    w =4 
elseif order == 2 
    xi = 1/sqrt(3).*s
    eta = 1/sqrt(3).*t 
    w = ones(1,4)
end
numIntegrationPoints = length(w)
#defining the material matrix plane subtriangulation
ν = 0.3
E = 210e9
ν′ = ν/(1-ν)
E′ = E/(1-ν^2)

D = E′/(1-ν′^2).*[1.0 ν′ 0; ν′ 1 0; 0 0 (1-ν′)/2]
#Compute the shape functions 
# B_mat = zeros(3,8,numIntegrationPoints,numElements)
# det_J_mat = zeros(numIntegrationPoints,numElements)
# for i_element 


