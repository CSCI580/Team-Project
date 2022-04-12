import random
from enum import Enum

class materialType(Enum):
    PHONGMATERIAL = "PHONGMATERIAL"
    BASICMATERIAL = "BASICMATERIAL"
    LAMBERTMATERIAL = "LAMBERTMATERIAL"
    
class reflectType(Enum):
    NONEREFLECT = "NONEREFLECT"
    MIRRORREFLECT = "MIRRORREFLECT"
    GLASSREFLECT = "GLASSREFLECT"

class textureType(Enum):
    NONE = "NONE"
    CHECKERBOARD = "CHECKERBOARD"
    MYSPECIAL = "MYSPECIAL"

class Material:
    data = dict()
        
    def __init__(self, color : list = [0.75,0.75,0.75], 
                 type : materialType = materialType.PHONGMATERIAL, 
                 shininess : float = 1000, 
                 specular : list = [1,1,1], 
                 reflectType: reflectType = reflectType.NONEREFLECT,
                 reflectivity : float = 1,
                 refractionRatio : float = 1.31,
                 texture : textureType = textureType.NONE) -> None:
       
        self.data = {
            "type": type,
            "color": color, #original color
            "shininess": shininess,  #determines how light difuses
            "specular": specular, #specular reflection color
            "reflectivity": reflectivity, #how reflective it is
            "reflectType": reflectType, 
            "refractionRatio": refractionRatio, #default 1.31 for ice
            "special": texture #texture
       }
      
    def __getitem__(self, key):
            return self.data[key]
    
    def __setitem__(self, key, val):
        self.data[key] = val
          
    def __iter__(self):
        for key in self.data:
            yield (key, self.data[key].value if isinstance(self.data[key],Enum) else self.data[key] )
                
class Object:
    data = dict()
    
    def __init__(self, t : str = "Object", mat : Material = Material(), com : str  = "Object" ) -> None:
        self.data = {
            "type" : t,
            "material" : mat,
            "comment" : "//"+com
        }
       
    def __getitem__(self, key):
        return self.data[key]
    
    def __setitem__(self, key, val):
        self.data[key] = val
     
    def __iter__(self):
        for key in self.data:
            if isinstance(self.data[key],Enum):
                yield (key, self.data[key].value )    
            elif isinstance(self.data[key],Material):
                yield (key, dict(self.data[key]) )    
            else:
                yield (key, self.data[key])    
        
class Plane(Object):
    def __init__(self, material: Material = Material(), normal: list = [0,0,1], 
                 dist: float = 0, comment: str = "Plane") -> None:
        super().__init__(t = "plane",mat = material, com = comment)
        self["normal"] = normal
        self["dist"] = dist
  
class Sphere(Object):
    def __init__(self, material: Material = Material(), center: list = [0,0,0], 
                 radius: float = 1, comment: str = "Sphere") -> None:
        super().__init__(t = "sphere", mat = material, com = comment)
        self["center"] = center
        self["radius"] = radius
  
class Box(Object):
    def __init__(self, material: Material = Material(), minCorner : list = [0,0,0], 
                 maxCorner : list = [1,1,1], comment: str = "Box") -> None:
        super().__init__(t = "box", mat = material, com = comment)
        self["minCorner"] = minCorner
        self["maxCorner"] = maxCorner
 
class Cylinder(Object):
    def __init__(self, material: Material = Material(), topCenter : list = [1,1,1], 
                 bottomCenter : list = [0,0,0], radius : float = 1,  comment: str = "Cylinder") -> None:
        super().__init__(t = "cylinder", mat = material, com = comment)
        self["topCenter"] = topCenter
        self["bottomCenter"] = bottomCenter
        self["radius"] = radius

class Cone(Object):
    def __init__(self, material: Material = Material(), topCenter : list = [1,1,1], 
                 bottomCenter : list = [0,0,0], radius : float = 1, comment: str = "Cone") -> None:
        super().__init__(t = "cone", mat = material, com = comment)
        self["topCenter"] = topCenter
        self["bottomCenter"] = bottomCenter
        self["radius"] = radius   
   
class Mesh(Object):
    def __init__(self, objfile : str, scale : float = 1, offset : list = [0,0,0], 
                 material: Material = Material(), commment: str = "Mesh"):
        super().__init__(t = "mesh", mat = material, com=commment)
        self["objfile"] = objfile
        self["scale"] = scale
        self["offset"] = offset

class Light:
    data = dict()
        
    def __init__(self, pos : list = [-20, 20, 10], color : list = [1,1,1],
                 intensity : float = 50, attenuate : float = 1) -> None:
        self.data = {
            "pos": pos,
            "color": color,
            "intensity": intensity, 
            "attenuate": attenuate,
        }
      
    def __getitem__(self, key):
        return self.data[key]
    
    def __setitem__(self, key, val):
        self.data[key] = val
          
    def __iter__(self):
        for key in self.data:
            yield (key, self.data[key])

