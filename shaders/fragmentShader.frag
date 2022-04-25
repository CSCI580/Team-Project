// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// flag for using soft shadows (set to 1 only when using soft shadows)
#define SOFT_SHADOWS 1

// define number of soft shadow samples to take
#define SOFT_SAMPLING 10

// define constant parameters
// EPS is for the precision issue
#define INFINITY 1.0e+12
#define EPS 1.0e-4
#define EPS_P 1.0e-8
#define M_PI 3.1415926535897932384626433832795

// define maximum recursion depth for rays
#define MAX_RECURSION 10

// define constants for scene setting
#define MAX_LIGHTS 10
#define AMBIENT vec3(0.2,0.2,0.2)

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3
#define ICE 4

struct Shape {
  int shapeType;
  vec3 v1;
  vec3 v2;
  float rad;
};

struct Material {
  int materialType;
  vec3 color;
  float shininess;
  vec3 specular;

  int materialReflectType;
  float reflectivity;
  float refractionRatio;
  int special;
};

struct Object {
  Shape shape;
  Material material;
};

struct Light {
  vec3 position;
  vec3 color;
  float intensity;
  float attenuate;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Intersection {
  vec3 position;
  vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

// varying
varying vec2 v_position;

// find then position some distance along a ray
vec3 rayGetOffset(Ray ray, float dist) {
  return ray.origin + (dist * ray.direction);
}


bool chooseCloserIntersection(float dist, inout float best_dist,
                              inout Intersection intersect,
                              inout Intersection best_intersect) {
  if (best_dist <= dist)
    return false;
  best_dist = dist;
  best_intersect.position = intersect.position;
  best_intersect.normal = intersect.normal;
  return true;
}


// Calculates surface normal of a triangle
vec3 findSurfaceNormal(vec3 t1, vec3 t2, vec3 t3) {
  vec3 u = t2 - t1;
  vec3 v = t3 - t1;

  vec3 normal = vec3(0.0, 0.0, 0.0);
  normal.xyz = u.yzx * v.zxy - u.zxy * v.yzx;

  return normalize(normal);
}

bool checkOneSide(vec3 v, vec3 P, vec3 t1, vec3 t2) {
  vec3 v1 = t1 - P;
  vec3 v2 = t2 - P;
  vec3 n1 = cross(v2, v1);

  return dot(v, n1) > EPS;
}

bool checkInBox(vec3 intersect, vec3 corner, vec3 size) {
  bool inside = true;
  intersect -= corner; // reframe coordinate
  for (int i = 0; i < 3; i++) {
    inside = inside && (intersect[i] > -EPS);
    inside = inside && (intersect[i] < size[i] + EPS);
  }
  return inside;
}


uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform float u_random;

#if 0
  float random (in vec2 st) {
      return fract(sin(dot(st.xy,
  vec2(12.9898,78.233)))*u_random);
  }
#else
  float random (in vec2 st) {
      return fract(sin(dot(st.xy,
  vec2(12.9898,78.233)))*
          43758.5453123);
  }
#endif

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise (in vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);

    // Four corners in 2D of a tile
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));

    vec2 u = f * f * (3.0 - 2.0 * f);

    return mix(a, b, u.x) +
            (c - a)* u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}

#define OCTAVES 6
float fbm (in vec2 st) {
    // Initial values
    float value = 0.0;
    float amplitude = .5;
    float frequency = 0.;
    //
    // Loop of octaves
    for (int i = 0; i < OCTAVES; i++) {
        value += amplitude * noise(st);
        st *= 2.;
        amplitude *= .5;
    }
    return value;
}


// forward declaration
float rayIntersectScene(Ray ray, out Material out_mat,
                        out Intersection out_intersect);

//Plane
float findIntersectionWithPlane(Ray ray, vec3 norm, float dist,
                                out Intersection intersect) {
  float a = dot(ray.direction, norm);
  float b = dot(ray.origin, norm) - dist;

  if (a < EPS && a > -EPS)
    return INFINITY;

  float len = -b / a;
  if (len < EPS)
    return INFINITY;

  intersect.position = rayGetOffset(ray, len);
  intersect.normal = norm;
  return len;
}

// Triangle
float findIntersectionWithTriangle(Ray ray, vec3 t1, vec3 t2, vec3 t3,
                                   out Intersection intersect) {

  vec3 faceNorm = findSurfaceNormal(t1, t2, t3);
  float d = dot(faceNorm, t1) / length(faceNorm); // proj
  float len = findIntersectionWithPlane(ray, faceNorm, d, intersect);

  if (!checkOneSide(ray.direction, intersect.position, t1, t2) ||
      !checkOneSide(ray.direction, intersect.position, t2, t3) ||
      !checkOneSide(ray.direction, intersect.position, t3, t1))
      return INFINITY;
  
  return len;

}

// Sphere
float findIntersectionWithSphere(Ray ray, vec3 center, float radius,
                                 out Intersection intersect) {

  vec3 L = center - ray.origin;
  float t_ca = dot(L, ray.direction);
  if (t_ca < EPS) return INFINITY;

  float d_square = dot(L, L) - t_ca * t_ca;
  if (d_square > radius * radius) return INFINITY;

  float t_hc = sqrt(radius * radius - d_square);

  float t1 = t_ca - t_hc;
  float t2 = t_ca + t_hc;

  if (t1 > EPS) {
    intersect.position = rayGetOffset(ray, t1);
    vec3 PO_diff = intersect.position - center;
    intersect.normal = PO_diff / length(PO_diff);
    return t1;
  }
  else if (t2 > EPS) {
    intersect.position = rayGetOffset(ray, t2);
    vec3 PO_diff = intersect.position - center;
    intersect.normal = PO_diff / length(PO_diff);
    return t2;
  }

  return INFINITY;

}

// Box
float findIntersectionWithBox(Ray ray, vec3 pmin, vec3 pmax,
                              out Intersection out_intersect) {

  vec3 size = pmax - pmin;
  vec3 norms[3];
  norms[0] = vec3(0.0, 0.0, 1.0); norms[1] = vec3(0.0, 1.0, 0.0); norms[2] = vec3(1.0, 0.0, 0.0);
  vec3 corners[2];
  corners[0] = pmin; corners[1] = pmax;

  Intersection intersect;
  float best_len = INFINITY;
  for (int i = 0; i < 6; i++) {
    vec3 norm = norms[int(mod(float(i), 3.0))];
    float dist = dot(norm, corners[int(i/3)]);
    float curr_len = findIntersectionWithPlane(ray, norm, dist, intersect);
    if (checkInBox(intersect.position, pmin, size)) {
      chooseCloserIntersection(curr_len, best_len, intersect, out_intersect);
    }
  }

  return best_len;

}

// Cylinder
float getIntersectOpenCylinder(Ray ray, vec3 center, vec3 axis, float len,
                               float rad, out Intersection intersect) {

  vec3 vd = ray.origin - center;
  float phi = dot(vd, axis);
  float theta = dot(ray.direction, axis);
  
  vec3 tmp1 = ray.direction - theta * axis;
  vec3 tmp2 = vd - phi * axis;

  float a = pow(length(tmp1), 2.0);
  float b = 2.0 * dot(tmp1, tmp2);
  float c = pow(length(tmp2), 2.0) - rad * rad;

  float t1 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
  float t2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);

  float res = INFINITY;

  if (!(t1 > -EPS && t2 > -EPS)) return INFINITY;
  if (t1 < t2 + EPS) {
    res = (t1 > -EPS) ? t1 : t2;
  } else{
    res = (t2 > -EPS) ? t2 : t1;
  }

  // find p2
  vec3 apex = center + axis * len;
  vec3 q = rayGetOffset(ray, res);
  if (dot(axis, q - center) > -EPS && dot(axis, q - apex) < EPS) {
    intersect.position = rayGetOffset(ray, res);
    vec3 PO_diff = intersect.position - center; 
    vec3 proj = dot(PO_diff, axis) * axis;
    vec3 normal = PO_diff - proj;
    intersect.normal = normal / length(normal);
    return res;
  }
  return INFINITY;

}

//disc
float getIntersectDisc(Ray ray, vec3 center, vec3 norm, float rad,
                       out Intersection intersect) {

  Intersection curr;
  float d = dot(norm, center) / length(norm);
  float len = findIntersectionWithPlane(ray, norm, d, curr);

  if (pow(length(curr.position - center), 2.0) < rad * rad) {
    intersect = curr;
    return len;
  }
  return INFINITY;

}

float findIntersectionWithCylinder(Ray ray, vec3 center, vec3 apex,
                                   float radius,
                                   out Intersection out_intersect) {
  vec3 axis = apex - center;
  float len = length(axis);
  axis = normalize(axis);

  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cylinder
  dist = getIntersectOpenCylinder(ray, center, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- two caps
  dist = getIntersectDisc(ray, center, -axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  dist = getIntersectDisc(ray, apex, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  return best_dist;
}

// Cone
float getIntersectOpenCone(Ray ray, vec3 apex, vec3 axis, float len,
                           float radius, out Intersection intersect) {

  vec3 vd = ray.origin - apex;
  float phi = dot(vd, axis);
  float theta = dot(ray.direction, axis);

  float side = sqrt(radius*radius + len*len);
  float cos_alpha2 = pow(len / side, 2.0);
  float sin_alpha2 = pow(radius / side, 2.0);

  vec3 tmp1 = ray.direction - theta * axis;
  vec3 tmp2 = vd - phi * axis;

  float a = pow(length(tmp1), 2.0) * cos_alpha2 - pow(theta, 2.0) * sin_alpha2;
  float b = 2.0 * (dot(tmp1, tmp2) * cos_alpha2 - theta * phi * sin_alpha2);
  float c = pow(length(tmp2), 2.0) * cos_alpha2 - pow(phi, 2.0) * sin_alpha2;

  float t1 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
  float t2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
  
  float res = INFINITY;

  if (!(t1 > -EPS && t2 > -EPS)) return INFINITY;
  if (t1 < t2 + EPS) {
    res = (t1 > -EPS) ? t1 : t2;
  } else{
    res = (t2 > -EPS) ? t2 : t1;
  }

  vec3 center = apex + axis * len;
  vec3 q = rayGetOffset(ray, res);
  if (dot(axis, q - apex) > -EPS && dot(axis, q - center) < EPS) {
    intersect.position = q;
    vec3 normal = cross(cross(axis, apex - q), apex - q);
    intersect.normal = normal/length(normal);
    return res;
  }

  return INFINITY;

}

float findIntersectionWithCone(Ray ray, vec3 center, vec3 apex, float radius,
                               out Intersection out_intersect) {
  vec3 axis = center - apex;
  float len = length(axis);
  axis = normalize(axis);

  Intersection intersect;
  float best_dist = INFINITY;
  float dist;


  dist = getIntersectOpenCone(ray, apex, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);


  dist = getIntersectDisc(ray, center, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  return best_dist;
}




bool pointInShadow(vec3 pos, vec3 lightVec) {

  Ray ray;
  ray.origin = pos;
  ray.direction = normalize(lightVec);
  

  Material hitMaterial;
  Intersection intersect;
  float len = rayIntersectScene(ray, hitMaterial, intersect);

  return (length(lightVec) - len > -EPS);

}


float softShadowRatio(vec3 pos, vec3 lightVec) {
  float not_in_shadow = 0.0;
  const int N = int(sqrt(float(SOFT_SAMPLING)));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // (i, j) represents tile locations
      // rand values are offsets in that grid, which has the range [0, 1/N]
      float x = float(i) / float(N) + random(vec2(0.0, 1.0/float(N)));
      float y = float(j) / float(N) + random(vec2(0.0, 1.0/float(N)));
      // convert to sphere coordinates
      float theta = 2.0 * M_PI * x;
      float phi = acos(2.0 * y - 1.0);
      float u = cos(phi);

      // find (x, y, z) offset
      vec3 sample = vec3(sqrt(1.0 - u*u)*cos(theta), sqrt(1.0 - u*u)*sin(theta), u);

      // apply sample offset to lightVec, check whether if the point is in shadow given the new light vec

      bool in_shadow = pointInShadow(pos, lightVec + sample);
      not_in_shadow += float(!in_shadow);
    }
  }
  return not_in_shadow/float(SOFT_SAMPLING);

}

vec3 getLightContribution(Light light, Material mat, vec3 posIntersection,
                          vec3 normalVector, vec3 eyeVector, bool phongOnly,
                          vec3 diffuseColor) {
  vec3 lightVector = light.position - posIntersection;

  if ((SOFT_SHADOWS == 0) && pointInShadow(posIntersection, lightVector)) return vec3(0.0, 0.0, 0.0);
  float ratio = SOFT_SHADOWS == 0 ? 1.0 : softShadowRatio(posIntersection, lightVector);


  if (ratio < EPS) return vec3(0.0, 0.0, 0.0);

  float distToLight = length(lightVector);
  lightVector /= length(lightVector);

  if (mat.materialType != PHONGMATERIAL && mat.materialType != LAMBERTMATERIAL)  return ratio * diffuseColor;

  vec3 contribution = vec3(0.0, 0.0, 0.0);

  // get light attenuation
  float attenuation = light.attenuate * distToLight;
  float diffuseIntensity =
      max(0.0, dot(normalVector, lightVector)) * light.intensity;

  // glass and mirror objects have specular highlights but no diffuse lighting
  if (!phongOnly) 
    contribution +=
        diffuseColor * diffuseIntensity * light.color / attenuation;
  
  if (mat.materialType == PHONGMATERIAL ) {
    vec3 n_eyeVector = normalize(eyeVector);
    vec3 R = -reflect(lightVector, normalVector);
    float intensity = pow(max(0.0, dot(n_eyeVector, R)), mat.shininess) * light.intensity;
    contribution += mat.specular * intensity * light.color / attenuation;
  }

  return ratio * contribution;

}



// find reflection or refraction direction (depending on material type)
vec3 calcReflectionVector(Material material, vec3 direction, vec3 normalVector,
                          bool isInsideObj) {

  if(material.materialReflectType == MIRRORREFLECT) 
    return reflect(direction, normalVector);
  
  vec3 random_vec = normalize(vec3(random(direction.xy + v_position), random(direction.yz + v_position), random(direction.zx + v_position)));
  random_vec = normalize(random_vec + normalVector);

  //non reflect 
  if (material.materialReflectType == NONEREFLECT ) 
    return normalize(reflect(direction, normalVector)+random_vec*0.35);


  float eta =
      (isInsideObj) ? (1.0 / material.refractionRatio) : material.refractionRatio;
  float cos_theta_i = dot(direction, normalVector)/(length(direction) * length(normalVector));
  float theta_i = acos(cos_theta_i);
  
  // reflect
  if ((eta * sin(theta_i)) > 1.0) 
    return normalize(reflect(direction, normalVector) + random_vec*0.05);
  
  // refract
  float theta_r = asin(eta * sin(theta_i));
  vec3 T = eta * direction - (eta * cos_theta_i + cos(theta_r)) * normalVector;

  return normalize(T+random_vec*0.05);
}


vec3 calculateDiffuseColor(Material mat, vec3 posIntersection,
                           vec3 normalVector) {

  if (mat.special == CHECKERBOARD) {

    float scale = 0.05; 
    vec3 color1 = vec3(0.3, 0.3, 0.3);
    vec3 color2 = vec3(0.5451, 0.2705, 0.0745);

    float a = floor(posIntersection.x * scale + EPS);
    float b = floor(posIntersection.y * scale + EPS);
    float c = floor(posIntersection.z * scale + EPS);

    vec3 res = (mod(a+b+c, 2.0) > 0.5) ? color1 : color2;
    return res;
  } 
  
  if (mat.special == MYSPECIAL) {

    float u_time = float(frame) * 0.1;
    float scale = 2.0;
    vec2 st = posIntersection.xy / scale;
    vec3 color = vec3(0.0);

    vec2 q = vec2(0.0);
    q.x = fbm(st + 0.0 * u_time);
    q.y = fbm(st + vec2(1.0));

    vec2 r = vec2(0.);
    r.x = fbm( st + 1.0*q + vec2(1.7,9.2)+ 0.15*u_time );
    r.y = fbm( st + 1.0*q + vec2(8.3,2.8)+ 0.126*u_time);

    float f = fbm(st+r);

    color = mix(vec3(0.101961,0.619608,0.666667),
                vec3(0.666667,0.666667,0.498039),
                clamp((f*f)*4.0,0.0,1.0));

    color = mix(color, vec3(0,0,0.164706), clamp(length(q),0.0,1.0));

    color = mix(color, vec3(0.666667,1,1), clamp(length(r.x),0.0,1.0));

    return vec3((f*f*f+.6*f*f+.5*f)*color);
  }

  return mat.color;

}



vec3 calculateColor(Material mat, vec3 posIntersection, vec3 normalVector,
                    vec3 eyeVector, bool phongOnly, vec3 nextColor, vec3 nectCoord) {
  
  vec3 diffuseColor = calculateDiffuseColor(mat, posIntersection, normalVector);

  vec3 outputColor = vec3(0.0, 0.0, 0.0);

  //ambient
  outputColor += AMBIENT * diffuseColor;

  for (int i = 0; i < MAX_LIGHTS; i++) {
    if (i >= numLights) break;
    vec3 contribution = getLightContribution(lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor);
    outputColor += contribution;
  }

  vec3 R = normalize(-reflect(nectCoord - posIntersection, normalVector));
  
  float intensity = pow(max(0.0, dot(eyeVector, R)), mat.shininess);
    
  outputColor += mat.specular * nextColor * intensity * 2.0 ;

  return outputColor;
}


vec3 traceRay(Ray ray) {
  vec3 resColor = vec3(0.0, 0.0, 0.0);
  float absorbtion = 0.1;
  vec3 resWeight = vec3(1.0, 1.0, 1.0);

  bool isInsideObj = false;


  Material mats[MAX_RECURSION];
  vec3 normals[MAX_RECURSION];
  bool isReflective[MAX_RECURSION];
  vec3 origin[MAX_RECURSION+1];
  vec3 resWeights[MAX_RECURSION];


  origin[0] = ray.origin;
  int hits = 0;

  for (int depth = 0; depth < MAX_RECURSION; depth++) {
    Material hitMaterial;
    Intersection intersect;
    float len = rayIntersectScene(ray, hitMaterial, intersect);
    if (abs(len) <= EPS || len >= INFINITY) break; 
   
    hits++;

    //store materials
    mats[depth] = hitMaterial;

    //store normals
    vec3 posIntersection = intersect.position;
    vec3 normalVector    = intersect.normal;
    vec3 eyeVector = normalize(ray.origin - posIntersection);
    if (dot(eyeVector, normalVector) < 0.0) {
        normalVector = -normalVector;
        isInsideObj = true;
    } else isInsideObj = false;
    
    normals[depth] = normalVector;

    //store reflective flags
    isReflective[depth] = (hitMaterial.materialReflectType == MIRRORREFLECT ||
                      hitMaterial.materialReflectType == GLASSREFLECT ||
                      hitMaterial.materialReflectType == ICE );
    //store intersections
    ray.origin = rayGetOffset(ray, len);
    origin[depth+1] = ray.origin;
    vec3 nextDir = calcReflectionVector(hitMaterial, ray.direction, normalVector, isInsideObj);
    ray.direction = normalize(nextDir);

    //store resWeights
    if(depth == 0)
      resWeights[depth] = vec3(1.0,1.0,1.0);
    else
      resWeights[depth] = hitMaterial.reflectivity * resWeights[depth-1]
        * (isReflective[depth-1] ? 0.9:0.1);

    
  } 

  // calculate color
  vec3 nextColor = vec3(0.0, 0.0, 0.0);
  for (int i = MAX_RECURSION; i>=0; i--) {
    if(i<=hits-2){
      vec3 currColor = calculateColor(mats[i], origin[i+1], normals[i], 
        normalize(origin[i] - origin[i+1]), isReflective[i], nextColor, origin[i+2]);
      nextColor = currColor * resWeights[i];
      resColor+=nextColor;
    }  
  }

  return resColor ;
}



void main() {
  float cameraFOV = 0.8;
  vec3 direction = vec3(v_position.x * cameraFOV * width / height,
                        v_position.y * cameraFOV, 1.0);
  Ray ray;
  ray.origin = vec3(uMVMatrix * vec4(camera, 1.0));
  ray.direction = normalize(vec3(uMVMatrix * vec4(direction, 0.0)));

  vec3 res = traceRay(ray);

  gl_FragColor = vec4(res.x, res.y, res.z, 1.0);
}
