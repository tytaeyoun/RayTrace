#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

Plano Convex Lens
Ray reaches FLAT surface first  | )
Shows :
    initial RMS value of a ray bundle
    z coordinate that gives the minimum RMS value (paraxial focus)
    minimum RMS value
    Plot of X-Y position of initial rays
    Plot of X-Y position or rays at the paraxial focus
    Plot of propagaion of rays in seen in X-Z plane   
'''
import matplotlib.pyplot as plt
import RayTrace as rt
'''
form Optical Element, flat surface facing the origin of the ray.
refractive index = 1.5168
z-intercepts = 100 , 105
curvature = -0.02
Aperature radius = 50
'''
flat = rt.SphericalRefraction(100, 0.0, 1, 1.5168, 50)
curved = rt.SphericalRefraction(105, -0.02, 1.5168, 1, 50)
Bundle = rt.RayBundle()

print("initial RMS of the bundle = ", Bundle.RMS())
'''
First plot : initial X-Y positions of the rays in he bundle
'''
fig = plt.figure(1)
for ray in Bundle.Bundl():
    #fig = plt.figure(1)
    x = ray.p()[0]
    y = ray.p()[1]
    z = ray.p()[2]
    plt.title('Ray Bundle in X-Y plane before propagation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.scatter(x, y, c = 'r')
    plt.grid(True)

'''
Propagate the bundle through the two surfaces
'''
Bundle.bundlepropagate([flat, curved])

'''
"for" loop used here to find the z coordinate of the paraxial focus.
approximate z value was found by observation of its propagation plot.
Variable "RMSmin" to always take the lower value of the two RMS values 
compared, and "zMin" taking the z-coordinate that gives the lower RMS value
'''
zMin = 0
RMSmin = 100
for i in range(240, 600, 1):
    Output1 = rt.OutputPlane(i)
    Bundle.bundlepropagate([Output1])
    if RMSmin > Bundle.RMS():
        RMSmin = Bundle.RMS()
        zMin = i

print("The minimum value of RMS is given at z = ", zMin)
print("Minimum RMS value is ", RMSmin)

'''
Bundle was propagated unil z = 600 by the for loop above.
minBundle stops propagating at the paraxial focus found above to allow 
plotting of X-Y location of the rays at the paraxial focus
'''
minBundle = rt.RayBundle()
MinOutput = rt.OutputPlane(zMin)
# MinOutput = rt.OutputPlane(401)
minBundle.bundlepropagate([flat, curved, MinOutput])

'''
Second plot : X-Y position of rays at the paraxial focus
'''
fig = plt.figure(2)
for ray in minBundle.Bundl():
    #fig = plt.figure(2)
    x = ray.p()[0]
    y = ray.p()[1]
    z = ray.p()[2]

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Ray Bundle in X-Y plane at the paraxial focus')
    plt.scatter(x, y, c = 'r')
    plt.grid(True)

'''
Third plot : X-Z plane showing the changes in the rays' positions and 
            directions as they propagate through different mediums
'''
# fig = plt.figure(3)
# for ray in Bundle.Bundl():
#     x = [positions[0] for positions in ray.vertices()]
#     y = [positions[1] for positions in ray.vertices()]
#     z = [positions[2] for positions in ray.vertices()]
#     plt.xlabel('x')
#     plt.ylabel('z')
#     plt.title('Bundle propagation. X-Z plane')
#     plt.plot(x, z, 'b-')
#     plt.grid(True)

plt.show()