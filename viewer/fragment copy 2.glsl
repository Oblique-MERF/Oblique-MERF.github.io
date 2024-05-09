vec2 rayAabbIntersection(vec3 aabbMin, vec3 aabbMax, vec3 origin,
                         vec3 invDirection) {
  vec3 t1 = (aabbMin - origin) * invDirection;
  vec3 t2 = (aabbMax - origin) * invDirection;
  vec3 tMin = min(t1, t2);
  vec3 tMax = max(t1, t2);
  return vec2(max(tMin.x, max(tMin.y, tMin.z)),
              min(tMax.x, min(tMax.y, tMax.z)));
}

vec3 brewerBlues(float value) {
    const vec3 blues[9] = vec3[](
        vec3(247.0/255.0, 251.0/255.0, 255.0/255.0),  // 最浅的蓝色
        vec3(222.0/255.0, 235.0/255.0, 247.0/255.0),
        vec3(198.0/255.0, 219.0/255.0, 239.0/255.0),
        vec3(158.0/255.0, 202.0/255.0, 225.0/255.0),
        vec3(107.0/255.0, 174.0/255.0, 214.0/255.0),
        vec3(66.0/255.0, 146.0/255.0, 198.0/255.0),
        vec3(33.0/255.0, 113.0/255.0, 181.0/255.0),
        vec3(8.0/255.0, 81.0/255.0, 156.0/255.0),
        vec3(8.0/255.0, 48.0/255.0, 107.0/255.0)   // 最深的蓝色
    );

    // 确保输入值在0.0到1.0的范围内
    value = clamp(value, 0.0, 1.0);

    // 计算插值索引
    float index = value * (9.f - 1.0);
    int i = int(index);
    int j = min(i + 1, 9 - 1);
    float t = index - float(i);

// 对颜色进行插值
    return mix(blues[i], blues[j], t);
}

#define SIGMOID(DTYPE) DTYPE sigmoid(DTYPE x) { return 1.0 / (1.0 + exp(-x)); }
SIGMOID(vec3)
SIGMOID(vec4)

#define DENORMALIZE(DTYPE)\
DTYPE denormalize(DTYPE x, float min, float max) {\
    return min + x * (max - min);\
}
DENORMALIZE(float)
DENORMALIZE(vec3)
DENORMALIZE(vec4)

float densityActivation(float x) { return exp(x - 1.0f); }

float densityToAlpha(float x, float stepSize) {
  return 1.0 - exp(-x * stepSize);
}

// Component-wise maximum
float max3 (vec3 v) {
  return max (max (v.x, v.y), v.z);
}

// Projective contraction
vec3 contract(vec3 x) {
  vec3 xAbs = abs(x);
  float xMax = max3(xAbs);
  if (xMax <= 1.0) {
    return x;
  }
  float scale = 1.0 / xMax;
  vec3 z = scale * x;
  // note that z.a = sign(z.a) where a is the the argmax component
  if (xAbs.x >= xAbs.y && xAbs.x >= xAbs.z) {
    z.x *= (2.0 - scale); // argmax = 0
  } else if (xAbs.y >= xAbs.x && xAbs.y >= xAbs.z) {
    z.y *= (2.0 - scale); // argmax = 1
  } else {
    z.z *= (2.0 - scale); // argmax = 2
  }
  return z;
}

// Inverse projective contraction
vec3 inverseContract(vec3 z) {
  vec3 zAbs = abs(z);
  float zMax = max3(zAbs);
  if (zMax <= 1.0) {
    return z;
  }
  float eps = 1e-6;
  float invZMax = max(eps, 2.0 - zMax);
  float scale = 1.0 / invZMax;
  vec3 x = scale * z;
  if (zAbs.x >= zAbs.y && zAbs.x >= zAbs.z) {
    x.x = sign(x.x) * scale; // argmax = 0
  } else if (zAbs.y >= zAbs.x && zAbs.y >= zAbs.z) {
    x.y = sign(x.y) * scale; // argmax = 1
  } else {
    x.z = sign(x.z) * scale; // argmax = 2
  }
  return x;
}

// Sorts an array of length 5 in-place. This is hardcoded to 5 since a ray
// traverses up to 5 quadrants.
void sort5(inout float[5] array, int arrayLength) {
  float t;
  for (int i = 0; i < arrayLength; ++i) {
    for (int j = i+1; j < arrayLength; ++j) {
      if (array[j] < array[i]) {
        t = array[i];
        array[i] = array[j];
        array[j] = t;
      }
    }
  }
}

// A solution is invalid if it does not lie on the plane or is outside of 
// the bounding box
#define INF 1e25
#define SOLUTION_CHECK(T, P, AXIS)\
q = contract(o + T.AXIS * d);\
if (abs(q.AXIS - P.AXIS) > eps || any(lessThan(q, aabbMin - eps)) ||\
    any(greaterThan(q, aabbMax + eps))) {\
  T.AXIS = -INF;\
}

// First checks wether the computed cancidate solutions are actually lieing on
// the bounding box. Then of all the valid intersections we return the one with
// the highest t-value (tMax).
// o: origin
// d: direction
// t0: candiate solutions for intersections with minimum YZ, XZ, XY planes
// t1: candiate solutions for intersections with maximum YZ, XZ, XY planes
// aabbMin: minimum of axis-aligned bounding bound
// aabbMin: maximum of axis-aligned bounding bound
float getTMax(vec3 o, vec3 d, vec3 t0, vec3 t1,
  vec3 aabbMin, vec3 aabbMax) {
  float eps = 1e-3;
  vec3 q;

  // Invalid solutions are set to -INF and therefore ignored.
  SOLUTION_CHECK(t0, aabbMin, x)
  SOLUTION_CHECK(t0, aabbMin, y)
  SOLUTION_CHECK(t0, aabbMin, z)
  SOLUTION_CHECK(t1, aabbMax, x)
  SOLUTION_CHECK(t1, aabbMax, y)
  SOLUTION_CHECK(t1, aabbMax, z)
  return max(max3(t0), max3(t1));
}

// The following functions compute intersections between rays and axis-aligned
// planes in contracted space.
// The seven functions correspond to seven cases assiociated with the seven
// quadrants present in projective contraction. The functions are derived 
// by solving contract(o+t*d) for t.
// o: origin
// d: direction
// p: x, y and z components define axis-aligned planes that the ray (o, d) is
//    intersected against
//    (x -> YZ-plane, y -> XZ-plane, z -> XY-plane)
vec3 h(vec3 o, vec3 d, vec3 p) {
  return (p - o) / d;
}

vec3 h0(vec3 o, vec3 d, vec3 p) {
  vec3 t;
  t.x = (1.0 / (2.0 - p.x) - o.x) / d.x;
  t.y = (o.y - p.y * o.x) / (p.y * d.x - d.y);
  t.z = (o.z - p.z * o.x) / (p.z * d.x - d.z);
  return t;
}

vec3 h1(vec3 o, vec3 d, vec3 p) {
  vec3 t;
  t.x = (o.x - p.x * o.y) / (p.x * d.y - d.x);
  t.y = (1.0 / (2.0 - p.y) - o.y) / d.y;
  t.z = (o.z - p.z * o.y) / (p.z * d.y - d.z);
  return t;
}

vec3 h2(vec3 o, vec3 d, vec3 p) {
  vec3 t;
  t.x = (o.x - p.x * o.z) / (p.x * d.z - d.x);
  t.y = (o.y - p.y * o.z) / (p.y * d.z - d.y);
  t.z = (1.0 / (2.0 - p.z) - o.z) / d.z;
  return t;
}

vec3 h3(vec3 o, vec3 d, vec3 p) {
  vec3 t;
  t.x = (1.0 / (-p.x - 2.0) - o.x) / d.x;
  t.y = -(o.x*p.y + o.y) / (d.x*p.y + d.y);
  t.z = -(o.x*p.z + o.z) / (d.x*p.z + d.z);
  return t;
}

vec3 h4(vec3 o, vec3 d, vec3 p) {
  vec3 t;
  t.x = -(o.y*p.x + o.x) / (d.y*p.x + d.x);
  t.y = (1.0 / (-p.y - 2.0) - o.y) / d.y;
  t.z = -(o.y*p.z + o.z) / (d.y*p.z + d.z);
  return t;
}

vec3 h5(vec3 o, vec3 d, vec3 p) {
  vec3 t;
  t.x = -(o.z*p.x + o.x) / (d.z*p.x + d.x);
  t.y = -(o.z*p.y + o.y) / (d.z*p.y + d.y);
  t.z = (1.0 / (-p.z - 2.0) - o.z) / d.z;
  return t;
}

// Intersects ray with all seven quadrants to obtain t-values at which the ray
// exits a quadrant. We need to know these t-values since whenever we 
// enter a new quadrant the origin and direction of the ray in contracted space
// needs to be recomputed.
float[5] findTraversedQuadrants(vec3 o, vec3 d, float near) {
  float[5] listQuadrantTMax = float[](INF, INF, INF, INF, INF); // Rays traverse up to 5 quadrants
  int numQuadrantsTraversed = 0;
  float c1 = 1.0 - 1e-5;
  float c2 = 2.0 - 1e-4;
  vec3 aabbMin;
  vec3 aabbMax;
  vec3 t0;
  vec3 t1;
  float tMax;

  // core region
  aabbMin = vec3(-1.0, -1.0, -1.0);
  aabbMax = vec3(1.0, 1.0, 1.0);
  t0 = h(o, d, aabbMin);
  t1 = h(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);

  // We discard intersections with quadrants that lie behind the camera
  // (tMax < near). When a quadrant is not traversed, getTMax returns -INF
  // and therefore this check also discards these values.
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  // argmax(|o+t*d|) = 0, o[0]+t*d[0] >= 0
  aabbMin = vec3( c1, -c1, -c1);
  aabbMax = vec3( c2,  c1,  c1);
  t0 = h0(o, d, aabbMin);
  t1 = h0(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  // argmax(|o+t*d|) = 1, o[1]+t*d[1] >= 0
  aabbMin = vec3(-c1, c1, -c1);
  aabbMax = vec3(c1, c2, c1);
  t0 = h1(o, d, aabbMin);
  t1 = h1(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  // argmax(|o+t*d|) = 2, o[2]+t*d[2] >= 0
  aabbMin = vec3(-c1, -c1, c1);
  aabbMax = vec3(c1, c1, c2);
  t0 = h2(o, d, aabbMin);
  t1 = h2(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  // argmax(|o+t*d|) = 0, o[0]+t*d[0] < 0
  aabbMin = vec3(-c2, -c1, -c1);
  aabbMax = vec3(-c1, c1, c1);
  t0 = h3(o, d, aabbMin);
  t1 = h3(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  // argmax(|o+t*d|) = 1, o[1]+t*d[1] < 0
  aabbMin = vec3(-c1, -c2, -c1);
  aabbMax = vec3(c1, -c1, c1);
  t0 = h4(o, d, aabbMin);
  t1 = h4(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  // argmax(|o+t*d|) = 2, o[2]+t*d[2] < 0
  aabbMin = vec3(-c1, -c1, -c2);
  aabbMax = vec3(c1, c1, -c1);
  t0 = h5(o, d, aabbMin);
  t1 = h5(o, d, aabbMax);
  tMax = getTMax(o, d, t0, t1, aabbMin, aabbMax);
  if (tMax >= near) {
    listQuadrantTMax[numQuadrantsTraversed] = tMax;
    numQuadrantsTraversed++;
  }

  sort5(listQuadrantTMax, numQuadrantsTraversed);
  return listQuadrantTMax;
}

struct BboxSetupResults {
  vec3 o; // ray origin in uncontracted space
  vec3 d; // ray direction in uncontracted space
  vec2 TMinMax; // uncontraction-space t-values at which the ray
  // enters or exits the bbox
};

// This function is called whenever we enter a bbox. We compute
// origin and direction of the ray in uncontracted space and compute for which
// t-value (in uncontracted space) the ray enters/exits the bbox
// tP and tQ are two world-space t-values that must lie within th entered
// bbox, i.e. (o+tP*d) and (o+tQ*d) must lie within the
// entered bbox.
BboxSetupResults bboxSetup(vec3 o, vec3 d, float tP, float tQ) {
  BboxSetupResults r;

  // Which quadrant did we enter?
  vec3 xP = o + tP * d;
  r.o = xP;
  r.d = d;
  vec3 aabbMin = vec3(-1.0, -1.0, -1.0);
  vec3 aabbMax = vec3(1.0, 1.0, 1.0);
  r.TMinMax = rayAabbIntersection(aabbMin, aabbMax, r.o, 1.0 / r.d);

  return r;
}

struct QuadrantSetupResults {
  vec3 oContracted; // ray origin in contracted space
  vec3 dContracted; // ray direction in contracted space
  vec2 quadrantTMinMaxContracted; // contraction-space t-values at which the ray
  // enters or exits the current quadrant
};

// This function is called whenever we enter a new quadrant. We compute
// origin and direction of the ray in contracted space and compute for which
// t-value (in contracted space) the ray enters/exits the quadrant
// tP and tQ are two world-space t-values that must lie within th entered
// quadrant, i.e. contract(o+tP*d) and  contract(o+tQ*d) must lie within the
// entered quadrant.
QuadrantSetupResults quadrantSetup(vec3 o, vec3 d, float tP, float tQ) {
  QuadrantSetupResults r;

  // Which quadrant did we enter?
  vec3 xP = o + tP * d;
  vec3 xAbs = abs(xP);
  float xMax = max3(xAbs);

  // Get the AABB of the quadrant the point x is in
  // Non-squash case, central quadrant:
  vec3 aabbMin = vec3(-1.0, -1.0, -1.0);
  vec3 aabbMax = vec3(1.0, 1.0, 1.0);
  if (xMax > 1.0) {
    // The point is inside in one of the outer quadrants ("squash zone")
    if (xAbs.x >= xAbs.y && xAbs.x >= xAbs.z) {
      aabbMin.x = xP.x > 0.0 ? 1.0 : -2.0; // argmax = 0
      aabbMax.x = xP.x > 0.0 ? 2.0 : -1.0;
    } else if (xAbs.y >= xAbs.x && xAbs.y >= xAbs.z) {
      aabbMin.y = xP.y > 0.0 ? 1.0 : -2.0; // argmax = 1
      aabbMax.y = xP.y > 0.0 ? 2.0 : -1.0;
    } else {
      aabbMin.z = xP.z > 0.0 ? 1.0 : -2.0; // argmax = 2
      aabbMax.z = xP.z > 0.0 ? 2.0 : -1.0;
    }
  }

  // Estimate the direction of the ray in contracted space by computing the
  // vector difference with two different t-values that are guanteed to
  // correspond to points within the current quadrant
  r.oContracted = contract(xP);
  vec3 zQ = contract(o + tQ * d);
  r.dContracted = normalize(zQ - r.oContracted);

  // When is the ray exiting the current quadrant? We need this value in
  // order to know when we enter a new quadrant or when to terminate ray marching.
  // Note that im findTraversedQuadrants word-space t-values are computed, while
  // we compute here contraction-space t-values. The world-space t-values are
  // needed to robustly obtain two points (tP and tQ) that are guranteed to lie
  // within a quadrant. With help of these values we can generate two points
  // in contracted space from which we can estimate the ray origin and direction
  // in contracted space. However, once we raymarch in contracted space we need
  // the contraction-space t-value to conveniently check whether we are still 
  // in the same quadrant. Alternatively, one could convert the contraction-
  // space point to a world-space point and estimate a world space t-value, but
  // this has been found to be numerically unstable.
  r.quadrantTMinMaxContracted =
      rayAabbIntersection(aabbMin, aabbMax, r.oContracted, 1.0 / r.dContracted);
  return r;
}


struct OccupancyQueryResults {
  bool inEmptySpace;
  float tBlockMax;
};

OccupancyQueryResults queryOccupancyGrid(
    vec3 z, vec3 minPosition, vec3 oContracted,
    vec3 invDContracted, highp sampler2DArray occupancyGrid,
    float voxelSizeOccupancy, vec2 gridSizeOccupancy,
    float voxelsize) {
  OccupancyQueryResults r;
  vec3 posOccupancy;
  vec3 blockMin;
  vec3 blockMax;
  vec2 occupancy;
  float occmin;
  float occmax;
  posOccupancy = (z - minPosition) / voxelSizeOccupancy;
  blockMin = floor(posOccupancy);
  blockMax = floor(posOccupancy) + 1.0;
  vec3[2] planeUv;
  planeUv[0] = vec3((blockMin.xz + blockMax.xz) * 0.5 / gridSizeOccupancy, 0.0);
  planeUv[1] = vec3((blockMin.xz + blockMax.xz) * 0.5 / gridSizeOccupancy, 1.0);
  occupancy.r = texture(
    occupancyGrid,
    planeUv[0]
  ).r;
  occupancy.g = texture(
    occupancyGrid,
    planeUv[1]
  ).r;
  occupancy.rg = (255.0 * occupancy.rg + zinit) * voxelsize + minPosition.y;
  occmax = occupancy.g;
  occmin = occupancy.r;
  if(occmin < z.y){
    r.inEmptySpace = (occmax < z.y);
    blockMin = blockMin * voxelSizeOccupancy + minPosition;
    blockMax = blockMax * voxelSizeOccupancy + minPosition;
    if(r.inEmptySpace){
      blockMin.y = occmax;
      blockMax.y = z.y;
    }
    else{
      blockMin.y = occmin;
      blockMax.y = z.y;
    }
    r.tBlockMax = rayAabbIntersection(blockMin, blockMax, oContracted, invDContracted).y;
  }
  else{
    r.inEmptySpace = true;
    // r.tBlockMax = 0.0;
    blockMin = blockMin * voxelSizeOccupancy + minPosition;
    blockMax = blockMax * voxelSizeOccupancy + minPosition;
    blockMin.y = z.y;
    blockMax.y = minPosition.y;
    r.tBlockMax = rayAabbIntersection(blockMin, blockMax, oContracted, invDContracted).y;
  }
  // r.inEmptySpace = false;
  // r.tBlockMax = 0.0;
  return r;
}


#define QUERY_OCCUPANCY_GRID(tBlockMax_L, occupancyGrid, voxelSizeOccupancy, gridSizeOccupancy)\
if (tunContracted > tBlockMax_L) {\
  occupancyQueryResults =\
    queryOccupancyGrid(z, minPosition, r.o, invDContracted,\
                        occupancyGrid, voxelSizeOccupancy, gridSizeOccupancy, occvoxelsize);\
  tBlockMax_L = occupancyQueryResults.tBlockMax;\
  if (occupancyQueryResults.inEmptySpace) {\
    if (occupancyQueryResults.tBlockMax == -1.0){\
      break;\
    }\
    tunContracted = max(tunContracted, tBlockMax_L) + 0.5 * stepSizeContracted;\
    continue;\
  }\
}

vec2 queryOccupancyGridmin(
    vec3 z, vec3 minPosition, vec3 oContracted,
    vec3 invDContracted, highp sampler2DArray occupancyGrid,
    float voxelSizeOccupancy, vec2 gridSizeOccupancy,
    float voxelsize) {
  OccupancyQueryResults r;
  vec3 posOccupancy;
  vec3 blockMin;
  vec3 blockMax;
  vec2 occupancy;
  float occmin;
  float occmax;
  posOccupancy = (z - minPosition) / voxelSizeOccupancy;
  blockMin = floor(posOccupancy);
  blockMax = floor(posOccupancy) + 1.0;
  vec3[2] planeUv;
  planeUv[0] = vec3((blockMin.xz + blockMax.xz) * 0.5 / gridSizeOccupancy, 0.0);
  planeUv[1] = vec3((blockMin.xz + blockMax.xz) * 0.5 / gridSizeOccupancy, 1.0);
  occupancy.r = texture(
    occupancyGrid,
    planeUv[0]
  ).r;
  occupancy.g = texture(
    occupancyGrid,
    planeUv[1]
  ).r;
  occupancy.rg = (255.0 * occupancy.rg + zinit) * voxelsize + minPosition.y;
  return occupancy.rg;
}


#define QUERY_OCCUPANCY_GRID_MIN(tBlockMax_L, occupancyGrid, voxelSizeOccupancy, gridSizeOccupancy)\
occupancyQueryResultsMin =\
  queryOccupancyGridmin(z, minPosition, r.o, invDContracted,\
                      occupancyGrid, voxelSizeOccupancy, gridSizeOccupancy, occvoxelsize);\


void main() {
  // See the DisplayMode enum at the top of this file.
  // Runs the full model with view dependence.
  const int DISPLAY_NORMAL = 0;
  // Disables the view-dependence network.
  const int DISPLAY_DIFFUSE = 1;
  // Only shows the latent features.
  const int DISPLAY_FEATURES = 2;
  // Only shows the view dependent component.
  const int DISPLAY_VIEW_DEPENDENT = 3;
  // Only shows the coarse block grid.
  const int DISPLAY_COARSE_GRID = 4;

  // Set up the ray parameters in world space..
  float nearWorld = nearPlane;
  vec3 originWorld = vOrigin;
  vec3 directionWorld = normalize(vDirection);

#ifdef USE_SPARSE_GRID
  ivec3 iGridSize = ivec3(round(sparseGridGridSize));
  int iBlockSize = int(round(dataBlockSize));
  ivec3 iBlockGridBlocks = (iGridSize + iBlockSize - 1) / iBlockSize;
  ivec3 iBlockGridSize = iBlockGridBlocks * iBlockSize;
  vec3 blockGridSize = vec3(iBlockGridSize);
#endif

  // float[5] listQuadrantTMax = findTraversedQuadrants(originWorld,
  //     directionWorld, nearWorld);
  
  float farWorld = rayAabbIntersection(minPosition, -1.0 * minPosition, originWorld, 1.0 / directionWorld).y;

  float tP = nearWorld;
  float tQ = mix(nearWorld, farWorld, 0.5);

  // QuadrantSetupResults r = quadrantSetup(originWorld, directionWorld, tP, tQ);
  BboxSetupResults r = bboxSetup(originWorld, directionWorld, tP, tQ);
  float tunContracted = 0.0;
  // int quadrantIndex = 1;

  float tBlockMax_L0 = -INF;
  float tBlockMax_L1 = -INF;
  float tBlockMax_L2 = -INF;
  float tBlockMax_L3 = -INF;
  float tBlockMax_L4 = -INF;

  float visibility = 1.0;
  vec3 accumulatedColor = vec3(0.0, 0.0, 0.0);
  vec4 accumulatedFeatures = vec4(0.0, 0.0, 0.0, 0.0);
  int step = 0;
  float astep = 0.f;
#ifdef USE_TRIPLANE
  #define GRID_SIZE triplaneGridSize
  #define VOXEL_SIZE triplaneVoxelSize
#else
  #define GRID_SIZE sparseGridGridSize
  #define VOXEL_SIZE sparseGridVoxelSize
#endif
  int maxStep = stepMult * int(ceil(length(GRID_SIZE)));
  float origStepSizeContracted = VOXEL_SIZE / float(stepMult);

  while (step < maxStep && visibility > 0.0 / 255.0) {
    step++;
#ifdef LARGER_STEPS_WHEN_OCCLUDED
    float stepSizeContracted = origStepSizeContracted *
        mix(8.0, 1.0, min(1.0, visibility / 0.66));
#else
    float stepSizeContracted = origStepSizeContracted;
#endif

    // check if the ray is exiting the current quadrant
    if (tunContracted > r.TMinMax.y) {
      break;
    }

    // Position of current sample in contracted space
    vec3 z = r.o + tunContracted * r.d;

    // Hierarchical empty space skipping
    vec3 invDContracted = 1.0 / r.d;
    OccupancyQueryResults occupancyQueryResults;
    QUERY_OCCUPANCY_GRID(tBlockMax_L0, occupancyGrid_L0, voxelSizeOccupancy_L0,
                         gridSizeOccupancy_L0)
    // QUERY_OCCUPANCY_GRID(tBlockMax_L1, occupancyGrid_L1, voxelSizeOccupancy_L1,
    //                      gridSizeOccupancy_L1)
    // QUERY_OCCUPANCY_GRID(tBlockMax_L2, occupancyGrid_L2, voxelSizeOccupancy_L2,
    //                      gridSizeOccupancy_L2)             
    // QUERY_OCCUPANCY_GRID(tBlockMax_L3, occupancyGrid_L3, voxelSizeOccupancy_L3,
    //                      gridSizeOccupancy_L3)
    QUERY_OCCUPANCY_GRID(tBlockMax_L4, occupancyGrid_L4, voxelSizeOccupancy_L4,
                         gridSizeOccupancy_L4)
 
    // We are in occupied space
    // compute grid positions for the sparse 3D grid and on the triplane planes
#ifdef USE_SPARSE_GRID
    vec3 posSparseGrid = (z - minPosition) / sparseGridVoxelSize - 0.5;
#endif
#ifdef USE_TRIPLANE
    vec3 posTriplaneGrid = (z - minPosition) / triplaneVoxelSize;
#endif

    // Calculate where the next sample would land in order to compute the
    // step size in world space (required for density-to-alpha conversion)
    // make sure not to shoot ouf the current quadrant
    float tunContractedNext = min(tunContracted + stepSizeContracted, r.TMinMax.y); 
    // Position of the next sample in contracted space
    vec3 zNext = r.o + tunContractedNext * r.d; 
    float stepSizeWorld = distance(zNext, z);

#ifdef USE_SPARSE_GRID
    vec3 atlasBlockMin =
        floor(posSparseGrid / dataBlockSize) * dataBlockSize;
    vec3 atlasBlockMax = atlasBlockMin + dataBlockSize;
    vec3 atlasBlockIndex =
        255.0 * texture(sparseGridBlockIndices, (atlasBlockMin + atlasBlockMax) /
                                      (2.0 * blockGridSize)).xyz;
    if (atlasBlockIndex.x <= 254.0) {
    vec3 posAtlas = clamp(posSparseGrid - atlasBlockMin, 0.0, dataBlockSize);

    posAtlas += atlasBlockIndex * (dataBlockSize + 1.0);
    posAtlas += 0.5;
    vec3 atlasUvw = posAtlas / atlasSize;

    if (displayMode == DISPLAY_COARSE_GRID) {
      // Half-pixel apron
      accumulatedColor = atlasBlockIndex * (dataBlockSize + 1.0) / atlasSize;
      // vec2 occupancyQueryResultsMin;
      // QUERY_OCCUPANCY_GRID_MIN(tBlockMax_L4, occupancyGrid_L4, voxelSizeOccupancy_L4,
      //                    gridSizeOccupancy_L4)

      // // float zzz = floor(occupancyQueryResultsMin.y*128.0)/128.0;
      // // if (z.y < occupancyQueryResultsMin.g && z.y > occupancyQueryResultsMin.r){
      // //   accumulatedColor = vec3(0.05,0.05,0.05);
      // // }
      // // float zzz = occupancyQueryResultsMin.r;
      // // accumulatedColor = vec3((zzz+1.0)/2.0,(zzz+1.0)/2.0,(zzz+1.0)/2.0);
      // float zzz = occupancyQueryResultsMin.g;
      // accumulatedColor = vec3((zzz+1.0)/2.0,(zzz+1.0)/2.0,(zzz+1.0)/2.0);
      accumulatedFeatures.rgb = atlasBlockIndex * (dataBlockSize + 1.0) / atlasSize;
      accumulatedFeatures.a = 1.0;
      visibility = 0.0;
      continue;
    }
#endif

    // First fetch all densities
#ifdef USE_SPARSE_GRID
    float density = texture(sparseGridDensity, atlasUvw).x;
    density = denormalize(density, rangeDensityMin, rangeDensityMax);
#else
    float density = 0.0;
#endif
#ifdef USE_TRIPLANE
    vec3[3] planeUv;
    planeUv[0] = vec3((posTriplaneGrid.y - triplaneResolutionZinit) / triplaneGridSizeZ, posTriplaneGrid.z / triplaneGridSize.x, 0.0);
    planeUv[1] = vec3(posTriplaneGrid.xz / triplaneGridSize, 0.0);
    planeUv[2] = vec3(posTriplaneGrid.x / triplaneGridSize.x, (posTriplaneGrid.y - triplaneResolutionZinit) / triplaneGridSizeZ, 0.0);

    float densityTemp;
    densityTemp = texture(planeDensity_1, planeUv[0]).x;
    densityTemp = denormalize(densityTemp, rangeDensityMin, rangeDensityMax);
    density += densityTemp;

    densityTemp = texture(planeDensity_2, planeUv[1]).x;
    densityTemp = denormalize(densityTemp, rangeDensityMin, rangeDensityMax);
    density += densityTemp;

    densityTemp = texture(planeDensity_3, planeUv[2]).x;
    densityTemp = denormalize(densityTemp, rangeDensityMin, rangeDensityMax);
    density += densityTemp;
#endif

    // Activate density and convert density to alpha.
    density = densityActivation(density);
    float alpha = densityToAlpha(density, stepSizeWorld);
    // if (z.y > -0.02)
    // { alpha = 0.0;}
    // Only fetch RGBFFFF (7 bytes) if alpha is non-negligible to save bandwidth
    if (alpha > 0.5 / 255.0) {
      astep = astep + 1.f;
#ifdef USE_SPARSE_GRID
      vec3 rgb = texture(sparseGridRgb, atlasUvw).rgb;
      rgb = denormalize(rgb, rangeFeaturesMin, rangeFeaturesMax);
#else
      vec3 rgb = vec3(0.0, 0.0, 0.0);
#endif
#ifdef USE_TRIPLANE
      vec3 rgbTemp;
      rgbTemp = texture(planeRgb_1, planeUv[0]).rgb;
      rgbTemp = denormalize(rgbTemp.rgb, rangeFeaturesMin, rangeFeaturesMax);
      rgb += rgbTemp;

      rgbTemp = texture(planeRgb_2, planeUv[1]).rgb;
      rgbTemp = denormalize(rgbTemp.rgb, rangeFeaturesMin, rangeFeaturesMax);
      rgb += rgbTemp;

      rgbTemp = texture(planeRgb_3, planeUv[2]).rgb;
      rgbTemp = denormalize(rgbTemp.rgb, rangeFeaturesMin, rangeFeaturesMax);
      rgb += rgbTemp;
#endif

      rgb = sigmoid(rgb); // Apply activation function

      if (displayMode != DISPLAY_DIFFUSE) {
        vec4 features = vec4(0.0, 0.0, 0.0, 0.0);
#ifdef USE_SPARSE_GRID
        features = texture(sparseGridFeatures, atlasUvw);
        features = denormalize(features, rangeFeaturesMin, rangeFeaturesMax);
#endif
#ifdef USE_TRIPLANE
        vec4 featuresTemp;
        featuresTemp = texture(planeFeatures_1, planeUv[0]);
        features +=
            denormalize(featuresTemp, rangeFeaturesMin, rangeFeaturesMax);

        featuresTemp = texture(planeFeatures_2, planeUv[1]);
        features +=
            denormalize(featuresTemp, rangeFeaturesMin, rangeFeaturesMax);

        featuresTemp = texture(planeFeatures_3, planeUv[2]);
        features +=
            denormalize(featuresTemp, rangeFeaturesMin, rangeFeaturesMax);
#endif

        features = sigmoid(features);
        accumulatedFeatures += visibility * alpha * features;
      }
      accumulatedColor += visibility * alpha * rgb;
      visibility *= 1.0 - alpha;
    }
#ifdef USE_SPARSE_GRID
    } // end of check: atlasBlockIndex.x <= 254.0
#endif
    tunContracted += stepSizeContracted;
  }

  if (displayMode == DISPLAY_VIEW_DEPENDENT) {
    accumulatedColor = vec3(0.0, 0.0, 0.0) * visibility;
  } else if (displayMode == DISPLAY_FEATURES) {
    accumulatedColor = accumulatedFeatures.rgb;
  }

  // Composite on white background
  accumulatedColor = vec3(1.0, 1.0, 1.0) * visibility + accumulatedColor;

  // Run view-dependency network
  if ((displayMode == DISPLAY_NORMAL ||
       displayMode == DISPLAY_VIEW_DEPENDENT)) {
    accumulatedColor += evaluateNetwork(accumulatedColor, accumulatedFeatures,
                             worldspaceROpengl * normalize(vDirection));
  }
  // accumulatedColor = brewerBlues(0.08*(astep-5.0)); 
  gl_FragColor = vec4(accumulatedColor, 1.0);
} 