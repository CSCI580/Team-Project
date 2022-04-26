import os
import json
import random
import math
from enum import Enum
from http.server import HTTPServer, CGIHTTPRequestHandler
from unittest import TextTestResult
from buildScene_util import *
if __name__ == "__main__":
    # generate objects
    objects = []
    
    # walls
    # p1 = Plane(comment = "back wall", normal=[0, 0, 1], dist = 100, 
    #            material=Material(color=[0.7, 0.7, 0.7],shininess=5,specular=[0.1,0.1,0.1]))
    # objects.append(p1)
    p2 = Plane(comment = "front wall", normal=[0, 0, -1], dist = 40, 
               material=Material(color=[0.7, 0.7, 0.7],shininess=5,specular=[0.1,0.1,0.1]))
    objects.append(p2)
    # p3 = Plane(comment = "right wall", normal=[1, 0, 0], dist = 40, 
    #            material=Material(color=[0.2, 0.7, 0.2],shininess=5,specular=[0.1,0.1,0.1]))
    # objects.append(p3)
    # p4 = Plane(comment = "left wall", normal=[-1, 0, 0], dist = 40, 
    #            material=Material(color=[0.7, 0.2, 0.2],shininess=5,specular=[0.1,0.1,0.1]))
    # objects.append(p4)
    # p5 = Plane(comment = "top wall", normal=[0, 1, 0], dist = 35, 
    #            material=Material(color=[0.7, 0.7, 0.7],shininess=5,specular=[0.1,0.1,0.1]))
    # objects.append(p5)
    # p6 = Plane(comment = "bottom wall", normal=[0, -1, 0], dist = 30, 
    #            material=Material(color=[0.7, 0.7, 0.7],shininess=50,texture=textureType.CHECKERBOARD))
    # objects.append(p6)

    

    # spheres
    
    s1 = Sphere(comment="bottom",center=[0,-835,30], radius=800, 
                material=Material(reflectType=reflectType.NONEREFLECT,color=[0.1, 0.1, 0.35], shininess=100, specular=[0.8,0.8,0.8]))
    objects.append(s1)
    s2 = Sphere(comment="left",center=[-840,0,30], radius=800, 
                material=Material(reflectType=reflectType.NONEREFLECT,color=[0.2, 0.7, 0.2], shininess=100, specular=[0.1,0.1,0.1]))
    objects.append(s2)
    s3 = Sphere(comment="right",center=[840,0,30], radius=800, 
                material=Material(reflectType=reflectType.NONEREFLECT,color=[0.7, 0.2, 0.2], shininess=100, specular=[0.1,0.1,0.1]))
    objects.append(s3)
    s4 = Sphere(comment="back",center=[0,0,900], radius=800, 
                material=Material(reflectType=reflectType.NONEREFLECT,color=[0.7, 0.7, 0.7], shininess=100, specular=[0.1,0.1,0.1]))
    objects.append(s4)
    s5 = Sphere(comment="top",center=[0,835,30], radius=800, 
                material=Material(reflectType=reflectType.NONEREFLECT,color=[0.2, 0.2, 0.2], shininess=100, specular=[0.1,0.1,0.1]))
    objects.append(s5)
    
    # s2 = Sphere(comment="the mirror sphere", center=[15, 2, 20], radius=8,
    #             material=Material(color=[0.2, 0.2, 0.2], shininess=1000, specular= [1, 1, 1], reflectivity=0.8, reflectType=reflectType.MIRRORREFLECT))
    # objects.append(s2)    
    
    
    # sg = Sphere(comment="the glass sphere", center=[0, 0, 30], radius=10,
    #             material=Material(color = [0.5,0.5,0.5],specular=[0.3,0.3,0.3], shininess=100, reflectivity=0.9, reflectType=reflectType.GLASSREFLECT, refractionRatio=1.309))
    # objects.append(sg)
    
    bubble_list = []
    with open("bubbleLocation.txt", 'r') as f:
        data = f.read().splitlines()
        for d in data:
            center = d.replace("[", "").replace(
                "]", "").split(', ')
            center = [float(i) for i in center]
            bubble_list.append(center)
    recentered_bubble_list = []
    random.shuffle(bubble_list)
    for i in range(50):
        vertex = bubble_list[i]
        vertex[0] = (vertex[0]/512*18 - 10) * 2.2 - 0.5
        vertex[1] = (vertex[1]/512*18 - 10) * 2.2 - 18
        vertex[2] = (vertex[2]/512*18 - 7) * 1.7 + 18

        recentered_bubble_list.append(vertex)
    # print(recentered_bubble_list)

    for center in recentered_bubble_list:
        s = Sphere(comment="random"+str(center),
                   center=center,
                   radius=random.uniform(0, 1)/3,
                   material=Material(color = [0,0,0],specular=[0,0,0], shininess=100, reflectivity=1, reflectType=reflectType.GLASSREFLECT, refractionRatio=1)
                   )
        objects.append(s)
    
    #cylinder
    # cy = Cylinder(comment="the cylinder", bottomCenter=[15.0, -30.0, 20.0], topCenter=[15.0, -6.0, 20.0], radius = 7,
    #               material=Material(shininess=100,texture=textureType.CHECKERBOARD))
    # objects.append(cy)
    
    # #cone
    # co = Cone(comment="the cone", bottomCenter=[-3, -5, 16], topCenter=[-3, 10, 16], radius=5,
    #           material=Material(color=[0.2, 0.4, 0.5],specular=[0.8, 0.8, 0.8], reflectType=reflectType.FRESNEL_BLEND_REFLECTION))
    # objects.append(co)
    
    # box
    # b = Box(comment="the box", minCorner= [-10, -30, 20], maxCorner=[10, -10, 40],
    #         material=Material(color = [0.5,0.5,0.5],specular=[0.3,0.3,0.3], shininess=100, reflectivity=0.9, reflectType=reflectType.GLASSREFLECT, refractionRatio=1.309))
    # objects.append(b)
    
    # b = Box(comment="the box", minCorner= [-10, -30, 20], maxCorner=[10, -10, 40],
    #         material=Material(reflectType=reflectType.ICE))
    # objects.append(b)
   
    drag = Mesh(objfile="./scenes/ice-sub-lowlow.obj", scale=12, offset=[-8,-25,30],
                material=Material(color = [0.5,0.5,0.5],specular=[0.3,0.3,0.3], shininess=100, reflectivity=0.9, reflectType=reflectType.GLASSREFLECT, refractionRatio=1.309))
    objects.append(drag)
    
    #generate lights
    lights = []
    l1 = Light(pos = [-30, 20, 60], color=[1, 1, 1], intensity=40, attenuate=5)
    lights.append(l1)
    l3 = Light(pos = [30, 20, 60], color=[1, 1, 1], intensity=40, attenuate=5)
    lights.append(l3)
    
    l2 = Light(pos = [0, 0, -10], color=[1, 1, 1], intensity=40, attenuate=1.5)
    lights.append(l2) 
   
    # l4 = Light(pos = [-30, -10, 30], color=[1, 1, 1], intensity=100, attenuate=1.5)
    # lights.append(l4) 
    
    #write to json file
    data = {"objects":[],"lights":[]}
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
