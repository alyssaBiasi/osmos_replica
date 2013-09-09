/***********************************************************************
 *      Realtime Rendering  -   Assignment 3
 *
 *       Alyssa Biasi s3328976
 *          and 
 *       Christopher Stojanovic s3334231
 *
 ***********************************************************************/

#ifndef UTILS_H
#define UTILS_H

#ifdef __APPLE__
#   include <OpenGL/gl.h>
#   include <OpenGL/glu.h>
#   include <GLUT/glut.h>
#else
#   include <GL/gl.h>
#   include <GL/glu.h>
#   include <GL/glut.h>
#endif

#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>

#define DEBUG
#undef DEBUG

/* Some handy macros */
#ifndef WIN32
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#define clamp(x, a, b) min(max(x, a), b)

#define PLAYER 0
#define AI 1
#define INIT_PARTICLE_SIZE 15
#define MIN_ZOOM 0.5

#define WIN 1
#define LOSE_SMALLEST -1
#define LOSE_ABSORBED 2
#define PAUSE 3
#define PLAY 0

#define BOUND 2.5
#define COOLTIME 1000
 
/* Debugging controls .*/
enum debugFlags {
  debug_time,
  debug_wall,
  debug_initialise_particle,
  debug_particle,
  debug_particle_collision,
  debug_collideParticlesBruteForce,
  debug_collideParticlesUniformGrid,
  debug_collisionReactionParticles2DbasisChange,
  debug_collisionReactionParticles2DprojNormal,
  debug_framerate,
  debug_range_check,
  debug_sum_kinetic_energy,
  debug_sum_momentum,
  numDebugFlags
};


/* Use our type so we can change precision easily. */
typedef double Real;

/* Particles (particles). */
struct Particle {
  Real position[2];
  Real velocity[2];
  Real radius;
  Real mass;
  Real elasticity;
  GLUquadric *quadric;  /* For rendering. */
  int slices, loops;    /* For rendering. */
  bool collided;
  int collidedWith;
  bool absorbed;
};


/* Control random or explicit initial positions */
enum {
  randomly,
  explicitly
} const initialiseParticles = randomly;


/* Collision detection method. */
enum CollisionDetectionMethod {
  bruteForce,
  uniformGrid
};


/* Control collision reaction calculation */
enum ReactionCalculation {
  basisChange,
  projNormal
} reacCalc = basisChange;


/* Arena. */
typedef struct {
  Real min[2], max[2];
  Real momentum[2];
  Real outline[360][2];
  int outcome;
  int level;
  Real dt;
  Real maxVelocity;
} Arena;


typedef struct {
  int width, height;
  float axisMax;
}Window;


/* The Camera struct is used for viewing the contents of the scene */
typedef struct {
  float fov;    /* Field of view of the camera, in degrees, in the y direction */
  float clipNear;   /* Depth value of the near clipping plane (not position on z axis!) */
  float clipFar;    /* Depth value of the far clipping plane (not position on z axis!) */
  bool zooming;   /* Flag to determine if camera is zooming */
  float zoom;   /* Distance camera is zoomed in/out */
  float scale;
  float maxZoom;
  float direction;
  float sensitivity;  /* Speed of camera rotation/zoom */
  Real pos[3];
} Camera;


#endif

