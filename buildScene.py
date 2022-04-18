import os
import json
import random
from enum import Enum
from http.server import HTTPServer, CGIHTTPRequestHandler
from buildScene_util import *
if __name__ == "__main__":
    # generate objects
    objects = []
    # walls
    p1 = Plane(comment="back wall", normal=[0, 0, 1], dist=40,
               material=Material(color=[0.2, 0.5, 0.2], texture=textureType.MYSPECIAL))
    objects.append(p1)
    p2 = Plane(comment="front wall", normal=[0, 0, -1], dist=40,
               material=Material(color=[0.5, 0.4, 0.1]))
    objects.append(p2)
    p3 = Plane(comment="right wall", normal=[1, 0, 0], dist=40,
               material=Material(color=[0.1, 0.1, 0.6]))
    objects.append(p3)
    p4 = Plane(comment="left wall", normal=[-1, 0, 0], dist=40,
               material=Material(color=[0.6, 0.1, 0.1]))
    objects.append(p4)
    p5 = Plane(comment="top wall", normal=[0, 1, 0], dist=35,
               material=Material(color=[0.3, 0.3, 0.3]))
    objects.append(p5)
    p6 = Plane(comment="bottom wall", normal=[0, -1, 0], dist=30,
               material=Material(color=[0.7, 0.7, 0.7], shininess=500, texture=textureType.CHECKERBOARD))
    objects.append(p6)
    # #spheres
    # s1 = Sphere(comment="the mat sphere",center=[-4, -22, 20], radius=8,
    #             material=Material(color=[0.2, 0.4, 0.5], shininess=200, specular=[0.8,0.8,0.8]))
    # objects.append(s1)
    # s2 = Sphere(comment="the mirror sphere", center=[15, 2, 20], radius=8,
    #             material=Material(color=[0.2, 0.2, 0.2], shininess=1000, specular= [1, 1, 1], reflectivity=0.8, reflectType=reflectType.MIRRORREFLECT))
    # objects.append(s2)
    # s3 = Sphere(comment="the glass sphere", center=[-20, -2, 15], radius=8,
    #             material=Material(specular=[1,1,1], shininess=70, reflectivity=0.8, reflectType=reflectType.GLASSREFLECT, refractionRatio=1.31))
    # objects.append(s3)
    # current box is from (-10, -10, 0) to (5, 5, 15)
    # bubble cloud 1: center(0, 0, 5) radius(~0.2)
    center1 = [0, 0, 5]
    radius1 = random.random()/2+0.1
    # for i in range(100):
    #     s = Sphere(comment="random"+str(i),
    #                center=center1,
    #                radius=radius1,
    #                material=Material(color=[random.random(), random.random(), random.random()],
    #                                  reflectType=reflectType.GLASSREFLECT))
    #     objects.append(s)
    # #cylinder
    # cy = Cylinder(comment="the cylinder", bottomCenter=[15.0, -30.0, 20.0], topCenter=[15.0, -6.0, 20.0], radius = 7,
    #               material=Material(shininess=100,texture=textureType.CHECKERBOARD))
    # objects.append(cy)
    # #cone
    # co = Cone(comment="the cone", bottomCenter=[-3, -5, 16], topCenter=[-3, 10, 16], radius=5,
    #           material=Material(color=[0.2, 0.4, 0.5],specular=[0.8, 0.8, 0.8], reflectType=reflectType.NONEREFLECT))
    # objects.append(co)
    # box
    b = Box(comment="the box", minCorner=[-10, -10, 0], maxCorner=[5, 5, 15],
            material=Material(specular=[1, 1, 1], shininess=70, reflectivity=0.8, reflectType=reflectType.GLASSREFLECT, refractionRatio=1.31))
    objects.append(b)
    # generate lights
    lights = []
    l1 = Light(pos=[-20, 20, 10], color=[1, 1, 1], intensity=50, attenuate=1.5)
    lights.append(l1)
    l2 = Light(pos=[20, 20, 10], color=[1, 1, 1], intensity=30, attenuate=2)
    lights.append(l2)
    l3 = Light(pos=[-10, 20, -20], color=[1, 1, 1], intensity=30, attenuate=1)
    lights.append(l3)
    l4 = Light(pos=[10, 20, -20], color=[1, 1, 1], intensity=20, attenuate=1)
    lights.append(l4)
    # write to json file
    data = {"objects": [], "lights": []}
    for i in objects:
        data['objects'].append(dict(i))
    for i in lights:
        data['lights'].append(dict(i))
    with open('./scenes/myscene.json', 'w') as outfile:
        json.dump(data, outfile, indent=4)
    # #hold server
    # # Make sure the server is created at current directory
    # os.chdir('.')
    # # Create server object listening the port 80
    # server_object = HTTPServer(server_address=('', 8000), RequestHandlerClass=CGIHTTPRequestHandler)
    # # Start the web server
    # server_object.serve_forever()
