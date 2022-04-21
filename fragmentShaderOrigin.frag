// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// flag for using soft shadows (set to 1 only when using soft shadows)
#define SOFT_SHADOWS 1

#define PI 3.14159265359
#define TWO_PI 6.28318530718
#define INV_PI 0.31830988618
#define POSITIVE_INF 1e10
#define NEGATIVE_INF -1e10
// define number of soft shadow samples to take
#define SOFT_SAMPLING 100

// define constant parameters
// EPS is for the precision issue
#define INFINITY 1.0e+12
#define EPS 1.0e-3
#define M_PI 3.1415926535897932384626433832795

// define maximum recursion depth for rays
#define MAX_RECURSION 50

// define constants for scene setting
#define MAX_LIGHTS 10

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
#define FRESNEL_BLEND_REFLECTION 5

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

mat3 inverse(mat3 m) {
  float a00 = m[0][0], a01 = m[0][1], a02 = m[0][2];
  float a10 = m[1][0], a11 = m[1][1], a12 = m[1][2];
  float a20 = m[2][0], a21 = m[2][1], a22 = m[2][2];
  float b01 = a22 * a11 - a12 * a21;
  float b11 = -a22 * a10 + a12 * a20;
  float b21 = a21 * a10 - a11 * a20;
  float det = a00 * b01 + a01 * b11 + a02 * b21;
  return mat3(b01, (-a22 * a01 + a02 * a21), (a12 * a01 - a02 * a11),
              b11, (a22 * a00 - a02 * a20), (-a12 * a00 + a02 * a10),
              b21, (-a21 * a00 + a01 * a20), (a11 * a00 - a01 * a10)) / det;
}

float random (in vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}
vec3 random3() {
  return vec3(random(v_position),random(v_position),random(v_position));
}
vec2 random2() {
  return random3().xy;
}
float random() {
  return random3().x;
}
float cosTheta(in vec3 w) {
  return w.z;
}
float cosTheta2(in vec3 w) {
  return w.z * w.z;
}
float sinTheta2(in vec3 w) {
  return 1.0 - cosTheta2(w);
}
float sinTheta(in vec3 w) {
  return sqrt(sinTheta2(w));
}
float tanTheta(in vec3 w) {
  return sinTheta(w) / cosTheta(w);
}
float tanTheta2(in vec3 w) {
  return sinTheta2(w) / cosTheta2(w);
}
float cosPhi(in vec3 w) {
  float s = sinTheta(w);
  return (s == 0.0) ? 1.0 : clamp(w.x / s, -1.0, 1.0);
}
float sinPhi(in vec3 w) {
  float s = sinTheta(w);
  return (s == 0.0) ? 0.0 : clamp(w.y / s, -1.0, 1.0);
}
float cosPhi2(in vec3 w) {
  float c = cosPhi(w);
  return c * c;
}
float sinPhi2(in vec3 w) {
  float s = sinPhi(w);
  return s * s;
}

// find then position some distance along a ray
vec3 rayGetOffset(Ray ray, float dist) {
  return ray.origin + (dist * ray.direction);
}

// if a newly found intersection is closer than the best found so far, record
// the new intersection and return true; otherwise leave the best as it was and
// return false.
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
float beckmannMicrofacetDistribution(in vec3 wh, in float alphaX, in float alphaY) {
  float tan2 = tanTheta2(wh);
  if (tan2 > POSITIVE_INF) return 0.0;
  float cos4 = pow(cosTheta(wh), 4.0);
  return exp(-tan2 * (cosPhi2(wh) / (alphaX * alphaX) + sinPhi2(wh) / (alphaY * alphaY))) /
      (PI * alphaX * alphaY * cos4);
}
float beckmannDistributionWhPdf(in vec3 wh, in float alphaX, in float alphaY) {
  return beckmannMicrofacetDistribution(wh, alphaX, alphaY) * abs(cosTheta(wh));
}
// put any general convenience functions you want up here
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 118 lines of code.

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
  if (dot(v, n1) < EPS) {
    return false;
  }
  return true;
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
mat3 orthonormal(in vec3 z) {
  vec3 w = normalize(z);
  vec3 v = normalize(cross(w, abs(w.x) > 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0)));
  vec3 u = normalize(cross(v, w));
  return mat3(u, v, w);
}
vec3 randomWiByCosine() {
  vec2 r = random2();
  float r1 = r.x;
  float r2 = r.y;
  float z = sqrt(1.0 - r2);
  float phi = TWO_PI * r1;
  float x = cos(phi) * sqrt(r2);
  float y = sin(phi) * sqrt(r2);
  return vec3(x, y, z);
}
vec3 randomWhByGGXDistribution(in vec3 wo, in float alphaX, in float alphaY) {
  float r1 = random();
  float phi = atan(alphaY / alphaX * tan(2.0 * PI * r1 + 0.5 * PI));
  phi += r1 > 0.5 ? PI : 0.0;
  float sinP = sin(phi);
  float cosP = cos(phi);
  float alpha2 = 1.0 / (cosP * cosP / (alphaX * alphaX) + sinP * sinP / (alphaY * alphaY));
  float r0 = random();
  float tan2Theta = alpha2 * r0 / (1.0 - r0);
  float cosT = 1.0 / sqrt(1.0 + tan2Theta);
  float sinT = sqrt(max(0.0, 1.0 - cosT * cosT));
  vec3 wh = vec3(sinT * cos(phi), sinT * sin(phi), cosT);
  wh *= wo.z * wh.z < 0.0 ? -1.0 : 1.0;
  return wh;
}

vec3 randomWiFresnelBlend(in vec3 wo) {
  if (random() < 0.5) {
    vec3 wh = randomWhByGGXDistribution(wo, 0.1, 0.1);
    return reflect(-wo, wh);
  } else {
    return randomWiByCosine();
  }
}

float randomFresnelBlendReflectionPdf(in vec3 wo, in vec3 wi) {
  vec3 wh = normalize(wo + wi);
  return 0.5 * (beckmannDistributionWhPdf(wh, 0.1, 0.1) / (4.0 * dot(wo, wh)))
             +  (abs(cosTheta(wi) * INV_PI));

}
// https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
highp float rand(vec2 co)
{
    highp float a = 12.9898;
    highp float b = 78.233;
    highp float c = 43758.5453;
    highp float dt= dot(co.xy ,vec2(a,b));
    highp float sn= mod(dt,3.14);
    return fract(sin(sn) * c);
}

// https://www.shadertoy.com/view/Xsl3Dl
// vec2 hash( vec2 p ) // replace this by something better
// {
// 	p = vec2( dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)) );
// 	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
// }

// float noise( in vec2 p )
// {
//     const float K1 = 0.366025404; // (sqrt(3)-1)/2;
//     const float K2 = 0.211324865; // (3-sqrt(3))/6;

// 	vec2  i = floor( p + (p.x+p.y)*K1 );
//     vec2  a = p - i + (i.x+i.y)*K2;
//     float m = step(a.y,a.x); 
//     vec2  o = vec2(m,1.0-m);
//     vec2  b = a - o + K2;
// 	vec2  c = a - 1.0 + 2.0*K2;
//     vec3  h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
// 	vec3  n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
//     return dot( n, vec3(70.0) );
// }

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

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
// float sinTheta(in vec3 w, in vec3 n) {
//   return sqrt(1 - cosTheta2(w, n));
// }
// float cosTheta(in vec3 w, in vec3 n) {
//   return dot(normalize(w), normalize(n));
// }

// float cosTheta2(in vec3 w, in vec3 n) {
//   return cosTheta(w,n)*cosTheta(w,n) ;
// }

// float tanTheta2(in vec3 w, in vec3 n) {
//   return (1 - cosTheta2(w, n) ) / cosTheta2(w, n);
// }
// float cosPhi(in vec3 w, in vec3 n) {
//   float s = sinTheta(w, n);
//   return (s == 0.0) ? 1.0 : clamp(w.x / s, -1.0, 1.0);
// }
// float cosPhi2(in vec3 w) {
//   float c = cosPhi(w);
//   return c * c;
// }

vec3 schlickFresnel(in vec3 f0, in float cosine) {
  return f0 + (vec3(1.0) - f0) * pow(1.0 - cosine, 5.0);
}
// float beckmannMicrofacetDistribution(in vec3 wh, in float alphaX, in float alphaY, in vec3 n) {
//   float tan2 = tanTheta2(wh, n);
//   if (tan2 > POSITIVE_INF) return 0.0;
//   float cos4 = pow(cosTheta(wh, n), 4.0);
//   return exp(-tan2 * (cosPhi2(wh) / (alphaX * alphaX) + sinPhi2(wh) / (alphaY * alphaY))) /
//       (PI * alphaX * alphaY * cos4);
// }


// vec3 fresnelBlendReflectionBRDF(in vec3 wo, in vec3 wi, in vec3 rd, in vec3 rs, in vec3 n) {
//   vec3 diffuse = rd * ((28.0) / (23.0 * PI)) * (vec3(1.0) - rs)
//                * (1.0 - pow(1.0 - 0.5 * abs(cosTheta(wi, n)), 5.0))
//                * (1.0 - pow(1.0 - 0.5 * abs(cosTheta(wo, n)), 5.0));
//   vec3 wh = normalize(wo + wi);
//   vec3 specular = beckmannMicrofacetDistribution(wh, 0.1, 0.1,n) * schlickFresnel(rs, dot(wi, wh))
//                 / (4.0 * abs(dot(wh, wi)) * max(cosTheta(wi,n), cosTheta(wo,n)));
//   return diffuse + specular;
// }
// ----------- STUDENT CODE END ------------

vec3 fresnelBlendReflectionBRDF(in vec3 wo, in vec3 wi, in vec3 rd, in vec3 rs) {
  vec3 diffuse = rd * ((28.0) / (23.0 * PI)) * (vec3(1.0) - rs)
               * (1.0 - pow(1.0 - 0.5 * abs(cosTheta(wi)), 5.0))
               * (1.0 - pow(1.0 - 0.5 * abs(cosTheta(wo)), 5.0));
  vec3 wh = normalize(wo + wi);
  vec3 specular = beckmannMicrofacetDistribution(wh, 0.1, 0.1) * schlickFresnel(rs, dot(wi, wh))
                / (4.0 * abs(dot(wh, wi)) * max(cosTheta(wi), cosTheta(wo)));
  return diffuse + specular;
}
// forward declaration
float rayIntersectScene(Ray ray, out Material out_mat,
                        out Intersection out_intersect);

// Plane
// this function can be used for plane, triangle, and box
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
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 28 lines of code.
  vec3 faceNorm = findSurfaceNormal(t1, t2, t3);
  float d = dot(faceNorm, t1) / length(faceNorm); // proj
  float len = findIntersectionWithPlane(ray, faceNorm, d, intersect);

  if (!checkOneSide(ray.direction, intersect.position, t1, t2) ||
      !checkOneSide(ray.direction, intersect.position, t2, t3) ||
      !checkOneSide(ray.direction, intersect.position, t3, t1))
      return INFINITY;
  
  return len;
  // ----------- STUDENT CODE END ------------
}

// Sphere
float findIntersectionWithSphere(Ray ray, vec3 center, float radius,
                                 out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 25 lines of code.
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
  // ----------- STUDENT CODE END ------------
}

// Box
float findIntersectionWithBox(Ray ray, vec3 pmin, vec3 pmax,
                              out Intersection out_intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // pmin and pmax represent two bounding points of the box
  // pmin stores [xmin, ymin, zmin] and pmax stores [xmax, ymax, zmax]
  // ----------- Our reference solution uses 44 lines of code.
  
  // size of the box
  vec3 size = pmax - pmin;
  vec3 norms[3];
  norms[0] = vec3(0.0, 0.0, 1.0); norms[1] = vec3(0.0, 1.0, 0.0); norms[2] = vec3(1.0, 0.0, 0.0);
  vec3 corners[2];
  corners[0] = pmin; corners[1] = pmax;

  // loop through all sides of the box
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
  // ----------- STUDENT CODE END ------------
}


// Cylinder
float getIntersectOpenCylinder(Ray ray, vec3 center, vec3 axis, float len,
                               float rad, out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 33 lines of code.
  // currently reports no intersection
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
  // ----------- STUDENT CODE END ------------
}

float getIntersectDisc(Ray ray, vec3 center, vec3 norm, float rad,
                       out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 18 lines of code.
  // currently reports no intersection
  Intersection curr;
  float d = dot(norm, center) / length(norm);
  float len = findIntersectionWithPlane(ray, norm, d, curr);

  if (pow(length(curr.position - center), 2.0) < rad * rad) {
    intersect = curr;
    return len;
  }
  return INFINITY;
  // ----------- STUDENT CODE END ------------
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
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 45 lines of code.
  // currently reports no intersection
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
  // ----------- STUDENT CODE END ------------
}

float findIntersectionWithCone(Ray ray, vec3 center, vec3 apex, float radius,
                               out Intersection out_intersect) {
  vec3 axis = center - apex;
  float len = length(axis);
  axis = normalize(axis);

  // -- infinite cone
  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cone
  dist = getIntersectOpenCone(ray, apex, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- caps
  dist = getIntersectDisc(ray, center, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  return best_dist;
}

vec3 calculateSpecialDiffuseColor(Material mat, vec3 posIntersection,
                                  vec3 normalVector) {
  // ----------- STUDENT CODE BEGIN ------------
  if (mat.special == CHECKERBOARD) {
    // ----------- Our reference solution uses 7 lines of code.
    float scale = 0.1; 
    vec3 color1 = vec3(1.0, 0.5, 0.5);
    vec3 color2 = vec3(0.6, 0.6, 1.0);

    float a = floor(posIntersection.x * scale + EPS);
    float b = floor(posIntersection.y * scale + EPS);
    float c = floor(posIntersection.z * scale + EPS);

    vec3 res = (mod(a+b+c, 2.0) > 0.5) ? color1 : color2;
    return res;
  } else if (mat.special == MYSPECIAL) {
    // ----------- Our reference solution uses 5 lines of code.
    // -----1
    // float scale = 7.0;
    // vec3 p = floor(posIntersection * scale + EPS);
    // float n = noise(p.xy);
    // float rand1 = rand(vec2(0.0, 1.0));
    // float rand2 = rand(vec2(0.0, 1.0));
    // vec3 col = vec3(rand1 * n, rand2 * n, n);
    // return col;

    // -----2
    // https://thebookofshaders.com/13/
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

    color = mix(color,
                vec3(0,0,0.164706),
                clamp(length(q),0.0,1.0));

    color = mix(color,
                vec3(0.666667,1,1),
                clamp(length(r.x),0.0,1.0));

    return vec3((f*f*f+.6*f*f+.5*f)*color);
  }

  // If not a special material, just return material color.
  return mat.color;
  // ----------- STUDENT CODE END ------------
}

vec3 calculateDiffuseColor(Material mat, vec3 posIntersection,
                           vec3 normalVector) {
  // Special colors
  if (mat.special != NONE) {
    return calculateSpecialDiffuseColor(mat, posIntersection, normalVector);
  }
  return vec3(mat.color);
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
bool pointInShadow(vec3 pos, vec3 lightVec) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 15 lines of code.
  Ray ray;
  ray.origin = pos;
  ray.direction = normalize(lightVec);

  Material hitMaterial;
  Intersection intersect;
  float len = rayIntersectScene(ray, hitMaterial, intersect);

  if (length(lightVec) - len > -EPS) return true;
  
  return false;
  // ----------- STUDENT CODE END ------------
}

// use random sampling to compute a ratio that represents the
// fractional contribution of the light to the position pos.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
float softShadowRatio(vec3 pos, vec3 lightVec) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 19 lines of code.
  float not_in_shadow = 0.0;
  const int N = int(sqrt(float(SOFT_SAMPLING)));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // (i, j) represents tile locations
      // rand values are offsets in that grid, which has the range [0, 1/N]
      float x = float(i) / float(N) + rand(vec2(0.0, 1.0/float(N)));
      float y = float(j) / float(N) + rand(vec2(0.0, 1.0/float(N)));
      // convert to sphere coordinates
      float theta = 2.0 * M_PI * x;
      float phi = acos(2.0 * y - 1.0);
      float u = cos(phi);

      // find (x, y, z) offset
      vec3 sample = vec3(sqrt(1.0 - u*u)*cos(theta), sqrt(1.0 - u*u)*sin(theta), u);

      // apply sample offset to lightVec, check whether if the point is in shadow given the new light vec
      bool in_shadow = pointInShadow(pos, lightVec + sample);
      not_in_shadow += in_shadow ? 0.0 : 1.0;
    }
  }
  return not_in_shadow/float(SOFT_SAMPLING);
  // ----------- STUDENT CODE END ------------
}

vec3 getLightContribution(Light light, Material mat, vec3 posIntersection,
                          vec3 normalVector, vec3 eyeVector, bool phongOnly,
                          vec3 diffuseColor) {
  vec3 lightVector = light.position - posIntersection;


  float ratio = 1.0; // default to 1.0 for hard shadows
  if (SOFT_SHADOWS == 1) {
    // if using soft shadows, call softShadowRatio to determine
    // fractional light contribution
    ratio = softShadowRatio(posIntersection, lightVector);
  }
  else {
    // check if point is in shadow with light vector
    if (pointInShadow(posIntersection, lightVector)) {
      return vec3(0.0, 0.0, 0.0);
    }
  }

  // Slight optimization for soft shadows
  if (ratio < EPS) {
    return vec3(0.0, 0.0, 0.0);
  }


  // normalize the light vector for the computations below
  float distToLight = length(lightVector);
  lightVector /= distToLight;

  if (mat.materialType == PHONGMATERIAL ||
      mat.materialType == LAMBERTMATERIAL) {
    vec3 contribution = vec3(0.0, 0.0, 0.0);

    // get light attenuation
    float attenuation = light.attenuate * distToLight;
    float diffuseIntensity =
        max(0.0, dot(normalVector, lightVector)) * light.intensity;

    // glass and mirror objects have specular highlights but no diffuse lighting
    if (!phongOnly) {
      contribution +=
          diffuseColor * diffuseIntensity * light.color / attenuation;
    }

    if (mat.materialType == PHONGMATERIAL) {
      // Start with just black by default (i.e. no Phong term contribution)
      vec3 phongTerm = vec3(0.0, 0.0, 0.0);
      // ----------- STUDENT CODE BEGIN ------------
      // ----------- Our reference solution uses 4 lines of code.
      vec3 n_eyeVector = normalize(eyeVector);
      vec3 R = -reflect(lightVector, normalVector);
      float intensity = pow(max(0.0, dot(n_eyeVector, R)), mat.shininess) * light.intensity;
      phongTerm = mat.specular * intensity * light.color / attenuation;
      // ----------- STUDENT CODE END ------------
      contribution += phongTerm;
    }

    return ratio * contribution;
  } else {
    return ratio * diffuseColor;
  }
}

vec3 calculateColor(Material mat, vec3 posIntersection, vec3 normalVector,
                    vec3 eyeVector, bool phongOnly) {
  // The diffuse color of the material at the point of intersection
  // Needed to compute the color when accounting for the lights in the scene
  vec3 diffuseColor = calculateDiffuseColor(mat, posIntersection, normalVector);

  // color defaults to black when there are no lights
  vec3 outputColor = vec3(0.0, 0.0, 0.0);

  // Loop over the MAX_LIGHTS different lights, taking care not to exceed
  // numLights (GLSL restriction), and accumulate each light's contribution
  // to the point of intersection in the scene.
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  // Return diffuseColor by default, so you can see something for now.
  for (int i = 0; i < MAX_LIGHTS; i++) {
    if (i >= numLights) break;
    vec3 contribution = getLightContribution(lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor);
    outputColor += contribution;
  }
  // return diffuseColor;
  return outputColor;
  // ----------- STUDENT CODE END ------------
}

// find reflection or refraction direction (depending on material type)
vec3 calcReflectionVector(Material material, vec3 direction, vec3 normalVector,
                          bool isInsideObj) {
  if (material.materialReflectType == MIRRORREFLECT) {
    return reflect(direction, normalVector);
  }
  // else if(material.materialReflectType == FRESNEL_BLEND_REFLECTION){
  //       mat3 ortho = orthonormal(normalVector);
  //       vec3 wo = inverse(ortho) * (-direction);
  //       vec3 wi = randomWiFresnelBlend(wo);
  //       float p = randomFresnelBlendReflectionPdf(wo, wi);
  //       if (p < 0.0) return vec3(0.0);
  //       vec3 dir = ortho * wi;
  //       vec3 brdf = fresnelBlendReflectionBRDF(wo, wi, info.surfaceInfo.color, info.surfaceInfo.anotherColor);
  //       c *= brdf * dot(normalVector, dir) / p;
  //       ray = Ray(info.pos, dir);
  // }
  // If it's not mirror, then it is a refractive material like glass.
  // Compute the refraction direction.
  // See lecture 13 slide (lighting) on Snell's law.
  // The eta below is eta_i/eta_r.
  // ----------- STUDENT CODE BEGIN ------------
  float eta =
      (isInsideObj) ? 1.0 / material.refractionRatio : material.refractionRatio;
  // ----------- Our reference solution uses 5 lines of code.
  // Return mirror direction by default, so you can see something for now.
  float cos_theta_i = dot(direction, normalVector)/(length(direction) * length(normalVector));
  float theta_i = acos(cos_theta_i);
  if ((eta * sin(theta_i)) > 1.) return reflect(direction, normalVector);
  float theta_r = asin(eta * sin(theta_i));
  vec3 T = eta * direction - (eta * cos_theta_i + cos(theta_r)) * normalVector;
  return T;
  // ----------- STUDENT CODE END ------------
}

vec3 traceRay(Ray ray) {
  // Accumulate the final color from tracing this ray into resColor.
  vec3 resColor = vec3(0.0, 0.0, 0.0);
  vec3 c =  vec3(1.0, 1.0, 1.0);

  // Accumulate a weight from tracing this ray through different materials
  // based on their BRDFs. Initially all 1.0s (i.e. scales the initial ray's
  // RGB color by 1.0 across all color channels). This captures the BRDFs
  // of the materials intersected by the ray's journey through the scene.
  vec3 resWeight = vec3(1.0, 1.0, 1.0);

  // Flag for whether the ray is currently inside of an object.
  bool isInsideObj = false;

  // Iteratively trace the ray through the scene up to MAX_RECURSION bounces.
  for (int depth = 0; depth < MAX_RECURSION; depth++) {
    // Fire the ray into the scene and find an intersection, if one exists.
    //
    // To do so, trace the ray using the rayIntersectScene function, which
    // also accepts a Material struct and an Intersection struct to store
    // information about the point of intersection. The function returns
    // a distance of how far the ray travelled before it intersected an object.
    //
    // Then, check whether or not the ray actually intersected with the scene.
    // A ray does not intersect the scene if it intersects at a distance
    // "equal to zero" or far beyond the bounds of the scene. If so, break
    // the loop and do not trace the ray any further.
    // (Hint: You should probably use EPS and INFINITY.)
    // ----------- STUDENT CODE BEGIN ------------
    Material hitMaterial;
    Intersection intersect;
    // ----------- Our reference solution uses 4 lines of code.
    float len = rayIntersectScene(ray, hitMaterial, intersect);
    if (abs(len) <= EPS || len >= INFINITY) break; 
    // ----------- STUDENT CODE END ------------

    // Compute the vector from the ray towards the intersection.
    vec3 posIntersection = intersect.position;
    vec3 normalVector    = intersect.normal;

    vec3 eyeVector = normalize(ray.origin - posIntersection);

    // Determine whether we are inside an object using the dot product
    // with the intersection's normal vector
    if (dot(eyeVector, normalVector) < 0.0) {
        normalVector = -normalVector;
        isInsideObj = true;
    } else {
        isInsideObj = false;
    }

    // Material is reflective if it is either mirror or glass in this assignment
    bool reflective = (hitMaterial.materialReflectType == MIRRORREFLECT ||
                       hitMaterial.materialReflectType == GLASSREFLECT || hitMaterial.materialReflectType == FRESNEL_BLEND_REFLECTION);

    // Compute the color at the intersection point based on its material
    // and the lighting in the scene
    vec3 outputColor = calculateColor(hitMaterial, posIntersection,
      normalVector, eyeVector, reflective);

    // A material has a reflection type (as seen above) and a reflectivity
    // attribute. A reflectivity "equal to zero" indicates that the material
    // is neither reflective nor refractive.

    // If a material is neither reflective nor refractive...
    // (1) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (2) Then break the for loop (i.e. do not trace the ray any further).
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 4 lines of code.
    if (!reflective) {

      outputColor = outputColor * resWeight;
      resColor += outputColor;
      break;
    }
    // ----------- STUDENT CODE END ------------

    // If the material is reflective or refractive...
    // (1) Use calcReflectionVector to compute the direction of the next
    //     bounce of this ray.
    // (2) Update the ray object with the next starting position and
    //     direction to prepare for the next bounce. You should modify the
    //     ray's origin and direction attributes. Be sure to normalize the
    //     direction vector.
    // (3) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (4) Update the current weight using the material's reflectivity
    //     so that it is the appropriate weight for the next ray's color.
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 8 lines of code.
    if (reflective) {
     
      ray.origin = rayGetOffset(ray, len);

      if(hitMaterial.materialReflectType == FRESNEL_BLEND_REFLECTION){
          mat3 ortho = orthonormal(normalVector);
          vec3 wo = inverse(ortho) * (-ray.direction);
          vec3 wi = randomWiFresnelBlend(wo);
          float p = randomFresnelBlendReflectionPdf(wo, wi);
          if (p < 0.0) return vec3(0.0);
          vec3 dir = ortho * wi;
          vec3 rd = vec3(30.0/255.0, 30.0/255.0, 240.0/255.0);
          vec3 rs =vec3(60.0/255.0, 60.0/255.0, 60.0/255.0);
          vec3 brdf = fresnelBlendReflectionBRDF(wo, wi, rd, rs);
          
          ray.direction = normalize(dir);

          c *= brdf * dot(normalVector, dir) / p;
          outputColor = outputColor * resWeight  * c * max(0.0,dot(normalVector,wi));
          resColor += outputColor;

          resWeight *= hitMaterial.reflectivity;
      }
      else{

          vec3 nextDir = calcReflectionVector(hitMaterial, ray.direction, normalVector, isInsideObj);
          ray.direction = normalize(nextDir);

          outputColor = outputColor * resWeight;
          resColor += outputColor;

          resWeight *= hitMaterial.reflectivity;
        }
      }
  }

  return resColor;
}
// const propsFresnelBlendReflection = {
//         diffuse: [30, 30, 240],
//         specular: [60, 60, 60],
//       }
//       const guiFresnelBlendReflection = gui.addFolder('Fresnel Blend Reflection');
//       addColorToGui(guiFresnelBlendReflection, propsFresnelBlendReflection, 'diffuse', function(v) {
//         glsl.setUniform('uFresnelBlendReflectionDiffuse', '3fv', [v[0] / 255.0, v[1] / 255.0, v[2] / 255.0]);
//       });
//       addColorToGui(guiFresnelBlendReflection, propsFresnelBlendReflection, 'specular', function(v) {
//         glsl.setUniform('uFresnelBlendReflectionSpecular', '3fv', [v[0] / 255.0, v[1] / 255.0, v[2] / 255.0]);
//       });
//       glsl.setUniform('uFresnelBlendReflectionDiffuse', '3fv',
//         [propsFresnelBlendReflection.diffuse[0] / 255.0, propsFresnelBlendReflection.diffuse[1] / 255.0, propsFresnelBlendReflection.diffuse[2] / 255.0]);
//         glsl.setUniform('uFresnelBlendReflectionSpecular', '3fv',
//           [propsFresnelBlendReflection.specular[0] / 255.0, propsFresnelBlendReflection.specular[1] / 255.0, propsFresnelBlendReflection.specular[2] / 255.0]);

void main() {
  float cameraFOV = 0.8;
  vec3 direction = vec3(v_position.x * cameraFOV * width / height,
                        v_position.y * cameraFOV, 1.0);

  Ray ray;
  ray.origin = vec3(uMVMatrix * vec4(camera, 1.0));
  ray.direction = normalize(vec3(uMVMatrix * vec4(direction, 0.0)));

  // trace the ray for this pixel
  vec3 res = traceRay(ray);

  // paint the resulting color into this pixel
  gl_FragColor = vec4(res.x, res.y, res.z, 1.0);
}
