#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:28:07 2018

@author: tyk15
"""
import numpy as np
import matplotlib.pyplot as plt

def normalise(a):
    normal = np.sqrt(np.dot(a,a))
    ahat = a/normal
    return ahat

class Ray:
    '''
    class modelling a Ray, records a ray with its position and direction in 3D 
    as arrays. Includes propagation of the ray (appending new position and 
    direction to the position and direction arrays) and plotting the ray in 
    X-Z plane
    '''
    def __init__(self, p = [0,0,0], k = [0,0,0]):
        self._position = [np.array(p)]
        self._direction = [np.array(k)/np.linalg.norm(k)]
        
        if len(self._position[-1]) != 3:
            raise Exception("Check position array size")
        if len(self._direction[-1]) != 3:
            raise Exception("Check direction array size")
        
    def __repr__(self):
        return "p=%r, k=%r" %( self._position[-1], self._direction[-1] )

    '''returns the last position of the ray'''
    def p(self): 
        return self._position[-1]
    
    '''returns the last direction of the ray'''
    def k(self): 
        return self._direction[-1]
    
    '''adds new position and direction to the ray'''
    def append(self, newp, newk): 
        if len(newp) != 3:
            raise Exception("Check position array size")
        if len(newk) != 3:
            raise Exception("Check direction array size")
            
        self._position.append(np.array(newp))
        self._direction.append(np.array(newk)/np.linalg.norm(newk))
        
    '''returns ray's position history'''    
    def vertices(self): 
        return self._position
    
    def propagate(self, l):
        self.append(self.p() + l * self.k(), self.k())
        

class OpticalElement:
    def propagate_ray(self,ray):
        "propagate a ray through the optical element"
        raise NotImplementedError()
        
class SphericalRefraction(OpticalElement):
    '''
    Class modelling a curved lens - includes: 
        function for calculating the position of its intercept with a ray 
        function giving the ray a new direction
        function that propagates the ray from the intercept to a new direction
    '''
    def __init__(self, z0 =105, curvature = 0.0, n1 = 1, n2 = 1.5, APradius = 100 ):
        '''
        z0 : z-intercept of the lenz
        curvature : 1/curvradius
        n1 : refractive index of medium ray started at
        n2 : refractive index of medium ray is entering
        APradius : aperature radius
        origin : position of origin the circular lens
        '''
        self._z0 = z0
        self._curvature =  curvature
        if self._curvature == 0:
            self._curvradius = 0
        else:
            self._curvradius = 1 / curvature
        self._n1 = n1
        self._n2 = n2
        self._APradius = APradius
        self._origin = np.array([0, 0, self._z0 + self._curvradius])
        
    def z0(self):
        return self._z0
    def curvradius(self):
        return self._curvradius
    def APradius(self):
        return self._APradius
    def origin(self):
        return self._origin
    def n1(self):
        return self._n1
    def n2(self):
        return self._n2
        
    def __repr__(self):
        return "z0 = %g\rcurvature = %g \rn1 = %g \rn2 = %g \rAperatureRadius = %g \rOrigin =%r"%\
         (self._z0, self._curvature, self._n1, self._n2, self._APradius, self._origin)
                        
    def intercept(self, ray):
        '''Returns the position where the ray and the Lens meets'''
        if self._curvature == 0:
            l = (self._z0 - ray.p()[2])/(ray.k()[2])
        else:
            #r : vector from Sphere centre to ray's starting point
            r = self._origin - ray.p()
            rMag = np.linalg.norm(r) #magnitude of r
            determinant = np.dot(r,ray.k())**2 - rMag**2 + self._curvradius**2

            #No Intercepts
            if determinant < 0:
                print ("No intercepts")
                return None        

            else:            
                lPlus  = -np.dot(r, ray.k()) + np.sqrt(determinant)
                lMinus = -np.dot(r, ray.k()) - np.sqrt(determinant)

                if self._curvature > 0:
                    l = np.abs(lPlus)
                if self._curvature < 0:
                    l = np.abs(lMinus)
        
        inters = ray.p() + l * ray.k()

        if np.sqrt(inters[0]**2 + inters[1]**2) > self._APradius:
            return None  
        else:
            return inters
    
    def newk(self, ray):
        '''returns new direction of the ray following Snell's Law'''
        
        '''Find the normal to the surface'''
        if self.intercept(ray).any() == None:
            return None
        if self._curvature == 0:
            n = np.array([0, 0, -1])
        elif self._curvature < 0:
            n = [0,0,self._curvradius] - self.intercept(ray)
        elif self._curvature > 0:
            n = -1 * ([0,0,self._curvradius + self._z0] - self.intercept(ray))        
        n = n/np.linalg.norm(n)
        
        '''Find the angle of incidence, theta1'''
        if self._curvature == 0:
            theta1 = np.pi
        else:  
            angleX = np.dot(n, ray.k())
            angleY = np.linalg.norm(n) * np.linalg.norm( ray.k() )
            theta1 = np.arccos(angleX/angleY)
            if theta1 > np.pi/2:
                theta1 = np.pi - theta1

        
        #total internal reflection
        nratio = self._n1/self._n2
        if np.sin(theta1) > nratio:
            print ("total internal reflection")
            return None
        
        '''Find the angle of refraction, theta2'''
        if theta1 == np.pi:
            theta2 = np.pi
        else:
            theta2 = np.arcsin( nratio * np.sin(theta1) )
            if theta2 > np.pi/2:
                theta2 = np.pi - theta2
        
        '''Find the new direction following the Snell's Law'''       
        newk1 = nratio * ray.k()
        factor1 = nratio * np.cos(theta1)
        factor2 = factor1 - np.cos(theta2)
        newk2 = factor2 * n
        
        newdirection = newk1 + newk2
        newdirection = newdirection/np.linalg.norm(newdirection)
        return newdirection
                    
    def PropagateRay(self, ray):  
        intercept = self.intercept(ray)
        
        if intercept.any() == None:
            print("No Intercept")
            return None
                
        else:
            newdirection = self.newk(ray)
            if self.newk(ray).any() == None:
                newdirection = np.array([0,0,0])
            ray.append(intercept, newdirection)
        
        
class OutputPlane(OpticalElement):
    '''
    the X-Y plane lying on some z axis, where rays are shown. 
    NO INTERACTION WITH RAY
    the intercepts are used to show the paths of rays
    '''
    def __init__ (self, z = 305):
        self.z = z
        
    def __repr__(self):
        return "%s(z = %r)" % ('OutputPlane', self.z)
    
    def intercept(self, ray):
        l = self.z - ray.p()[2]
        factor =l * ray.k()
        inters = ray.p() + factor
        return inters
        
    def PropagateRay(self, ray):
        if self.intercept(ray).all() == None:
            print ("No intercept with the output plane")
            return None
        
        ray.append(self.intercept(ray), ray.k())
        
class RayBundle:
    '''
    Makes a Bundle of rays in circular positions, with all rays parallel to 
    z-axis. Parameter R is the distance of the rays from the z-axis.
    Includes :
        function to plot the bundle 
        function to propagate each rays through OpticalElements
        function to calculate the Root mean square of the rays
    '''
    def __init__(self):
        self.Bundle = []
        self.direction = [0,0,1]
        self.R = [0.0, 2, 4, 6, 8, 10, 12]
        for i in range(len(self.R)):
            for j in range(i * 10):
                t = j*(2 * np.pi / (i * 10))
                self.Bundle.append(Ray([self.R[i] * np.cos(t), self.R[i] * np.sin(t), 0], self.direction))
    
    def Bundl(self):
        return self.Bundle
    
    def bundlepropagate(self, OpticalElements):
        for element in OpticalElements:
            for ray in self.Bundl():
                element.PropagateRay(ray)

    def RMS(self):
        """ 
        calculates the RMS value of the rays from the z axis
        Used to estimate the size of the geometrical focus
        """
        N = len(self.Bundl())
        x = []
        y = []
        sqrtSums = []
        
        for rays in self.Bundl():
            x.append(rays.p()[0])
            y.append(rays.p()[1])
        
        for i in range(len(x)):
            sqrtSums.append(np.sqrt(x[i]**2 + y[i]**2))
        Total = np.sum(sqrtSums)
        rms = Total/N
        return rms