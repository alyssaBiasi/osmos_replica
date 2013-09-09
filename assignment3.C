/***********************************************************************
 *      Realtime Rendering  -   Assignment 3
 *
 *       Alyssa Biasi s3328976
 *          and 
 *       Christopher Stojanovic s3334231
 *
 ***********************************************************************/

/*
 * $Id: particles_2D.C,v 1.13 2012/09/24 00:57:22 gl Exp gl $
 */

#include "utils.h"
 
int debug[numDebugFlags] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/* Small number to handle numerical imprecision. */
const Real epsilon = 1.0e-6;

/* Control 1D or 2D */ 
const int dimension = 2;

/* Rendering info. */
enum renderMode { wire, solid };
static renderMode renMode = solid;
static Real elapsedTime = 0.0, startTime = 0.0;
static const int milli = 1000;
static bool go = false;

bool cooling = false;
int coolPeriod = 0;

Particle *particle = new Particle[INIT_PARTICLE_SIZE];
int particle_size = INIT_PARTICLE_SIZE;
int numParticles = 0;

CollisionDetectionMethod CDmethod = uniformGrid;
Camera camera;
Window window;
Arena arena;

Real random_uniform() {
  return rand()/(float)RAND_MAX;
}


void panic(const char *m) {
    printf("%s", m);
    exit(1);
}


void quit() {
    free(particle);
    exit(EXIT_SUCCESS);
}


void initialiseArena() {
#ifdef DEBUG
    printf("init arena\n");
#endif

    const Real radius = 20;
    double pi = 4.0 * atan(1.0);
    int angle = 1;
    int index = 0;

    go = false;
    
    arena.min[0] = -radius;
    arena.min[1] = -radius;
    arena.max[0] = radius;
    arena.max[1] = radius;

    arena.momentum[0] = 0.0;
    arena.momentum[1] = 0.0;

    arena.dt = 0.0;
    arena.maxVelocity = 1.6;

    /* Storing circle coordinates. */
    while(angle <= 360) {
        arena.outline[index][0] = radius*cosf(angle*pi/180);
        arena.outline[index][1] = radius*sinf(angle*pi/180);
        angle++;
        index++;
    }
}


/*******************************************************************
*
*   Particle functions
*
*******************************************************************/

/* Determines whether the player is biggest/smallest mote */
int calculateMoteSize() {
    int i;
    float biggestMote = particle[PLAYER+1].mass;
    float smallestMote = particle[PLAYER+1].mass;
    float totalRemaingMass = 0;
    
    /* Determines the mass of largest & smallest mote */
    for(i=PLAYER+1; i<numParticles; i++) {
        if(particle[i].mass > biggestMote) {
            biggestMote = particle[i].mass;
        }
        else if(particle[i].mass < smallestMote) {
            smallestMote = particle[i].mass;
        }
    }

    /* Compare to player's mass */
    if(particle[PLAYER].mass > biggestMote) {
        return WIN;
    }
    else if(particle[PLAYER].mass < smallestMote) {
        return LOSE_SMALLEST;
    }

    /* Determines if possible for player to win */
    for(i=PLAYER+1; i<numParticles; i++) {
        if(particle[PLAYER].mass >= particle[i].mass) {
            totalRemaingMass += particle[i].mass;
        }
    }

    if(totalRemaingMass + particle[PLAYER].mass <= biggestMote) {
        return LOSE_SMALLEST;
    }

    return 0;
}


void resizeParticleArray() {
#ifdef DEBUG
    printf("resize particle array \n");
#endif

    int i;

    /* Resize array */
    particle_size *= 2;
    particle = (Particle*)realloc(particle, sizeof(Particle) * particle_size);

    if(particle == NULL) {
        printf("Error - my init: initialising particle array.\n");
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    printf("end resize particle array \n");
#endif
}


void removeParticle(int i) { 
    /* Swaps AI particle with the one that absorbed it */
    if(arena.level == 2 && i == AI) {
        int collidedWith = particle[i].collidedWith;

        if(collidedWith != PLAYER) {
            particle[i] = particle[collidedWith];
            i = collidedWith;
        }
    }

    for(; i<numParticles-1; i++) {
        particle[i] = particle[i+1];
    }
    
    numParticles--;
}


void initPreventCollisions(int i) {
#ifdef DEBUG
    printf("prevent collisions\n");
#endif
    Real n[2], n_mag_sq, sum_radii, sum_radii_sq;
    Real radius, angle;
    double pi = 4.0*atan(1.0);
    bool collision = false;
    bool done = false;
    int j;

    while(!done) {
        radius = random_uniform() * (arena.max[0] - 1.3 * particle[i].radius);
        angle = random_uniform() * 360;

        particle[i].position[0] = radius*cosf(angle*pi/180);
        particle[i].position[1] = radius*sinf(angle*pi/180);

        /* Check for collision with existing particles. */
        collision = false;
        j = 0;
        while (!collision && j < i) {
            sum_radii = particle[i].radius + particle[j].radius + 2.0;
            sum_radii_sq = sum_radii * sum_radii;
            n[0] = particle[j].position[0] - particle[i].position[0];
            n[1] = particle[j].position[1] - particle[i].position[1];
            n_mag_sq = n[0] * n[0] + n[1] * n[1];
            
            if (n_mag_sq < sum_radii_sq)
                collision = true;
            else
                j++;
        }

        if (!collision) {
            done = true;
        }
    }
}


void addParticle(bool checkCollisions) {
    GLUquadric *quadric = gluNewQuadric();
    int i;

    /* Resizes array if necessary */
    if(numParticles >= particle_size) {
        resizeParticleArray();
    }

    i = numParticles;
    
    particle[i].velocity[0] = (random_uniform() - 0.5) * arena.maxVelocity - 0.6;
    particle[i].velocity[1] = (random_uniform() - 0.5) * arena.maxVelocity - 0.6;
    particle[i].mass = random_uniform() * 1.8;
    particle[i].radius = sqrt(particle[i].mass);
    particle[i].elasticity = 1.0;
    particle[i].quadric = quadric;
    particle[i].slices = 10;
    particle[i].loops = 3;
    particle[i].collided = false;
    particle[i].absorbed = false;
    particle[i].collidedWith = -1;

    if(checkCollisions) {
        initPreventCollisions(i);
    }

    if (debug[debug_initialise_particle]) {
        printf ("initialiseParticles: x %f y %f\n", 
                particle[i].position[0], particle[i].position[1]);
    }

    numParticles++;
}


void initialiseParticlesRandomly() {
    GLUquadric *quadric = gluNewQuadric();
    double pi = 4.0*atan(1.0);
    int i;

#ifdef DEBUG
    printf("random particle init\n");
#endif

    /* Initialise the player particle */
    particle[PLAYER].velocity[0] = 0.00;
    particle[PLAYER].velocity[1] = 0.00;
    particle[PLAYER].mass = 0.5;
    particle[PLAYER].radius = sqrt(particle[PLAYER].mass);
    particle[PLAYER].elasticity = 1.0;
    particle[PLAYER].quadric = quadric;
    particle[PLAYER].slices = 10;
    particle[PLAYER].loops = 3;
    particle[PLAYER].collided = false;
    particle[PLAYER].absorbed = false;
    particle[PLAYER].collidedWith = -1;

    /* Random radius & angle to find position in circular arena. */
    random_uniform();
    Real radius = random_uniform() * (arena.max[0] - 2.0 * particle[PLAYER].radius);
    Real angle = random_uniform() * 360;

    particle[PLAYER].position[0] = radius*cosf(angle*pi/180);
    particle[PLAYER].position[1] = radius*sinf(angle*pi/180);
    numParticles++;

    for (i = PLAYER+1; i < INIT_PARTICLE_SIZE; i++) {
        /* Init all other particles to random values */
        addParticle(true);
    }

    /* Make sure that player isn't biggest/smallest mote. */
    while(calculateMoteSize() != 0) {
#ifdef DEBUG
    printf("resetting field\n");
#endif
        for(i=PLAYER+1; i < numParticles; i++) {
            particle[i].mass = random_uniform() * 1.0;
            particle[i].radius = sqrt(particle[i].mass);
        }
    }

    for(i=1; i<numParticles; i++) {
        initPreventCollisions(i);
    }

#ifdef DEBUG
    printf("end particle random init\n");
#endif
}


void propelPlayer(Real objX, Real objY) {
    int index = numParticles;
    float acceleration = 1.2;
    float dampening = 0.3;
    Real vecX, vecY;

    /* Adds ejected mass to array */
    addParticle(false);
    particle[index].mass = 0.0005;
    particle[index].radius = sqrt(particle[index].mass);
    particle[index].position[0] =  particle[PLAYER].position[0];
    particle[index].position[1] =  particle[PLAYER].position[1];

    /* Vector from player to click */
    vecX = objX - particle[PLAYER].position[0];
    vecY = objY - particle[PLAYER].position[1];

    particle[index].velocity[0] = vecX * acceleration;
    particle[index].velocity[1] = vecY * acceleration;

    /* Moves mass along path - no collision with player */
    particle[index].position[0] += particle[index].velocity[0] * 0.6;
    particle[index].position[1] += particle[index].velocity[1] * 0.6;

    /* Applying propellsion to player particle */
    particle[PLAYER].mass -= particle[index].mass;
    particle[PLAYER].radius = sqrt(particle[PLAYER].mass);
    particle[PLAYER].velocity[0] -= particle[index].velocity[0] * particle[PLAYER].mass * dampening;
    particle[PLAYER].velocity[1] -= particle[index].velocity[1] * particle[PLAYER].mass * dampening;
}


void propelAIAway(Real objX, Real objY) {
    int index = numParticles;
    float acceleration = 1.2;
    float dampening = 0.5;
    Real vecX, vecY;
    
    /* "Chaos" value to aid in escape */
    float chaos = (float)((rand() % 400) / 100) - 2.0;

    /* Adds ejected mass to array */
    addParticle(false);
    particle[index].mass = 0.0005;
    particle[index].radius = sqrt(particle[index].mass);
    particle[index].position[0] =  particle[AI].position[0];
    particle[index].position[1] =  particle[AI].position[1];

    /* Vector from AI to click */
    vecX = (objX + chaos) - particle[AI].position[0];
    vecY = (objY + chaos) - particle[AI].position[1];

    particle[index].velocity[0] = vecX * acceleration;
    particle[index].velocity[1] = vecY * acceleration;

    /* Moves mass along path - no collision with player */
    particle[index].position[0] += particle[index].velocity[0] * 0.6;
    particle[index].position[1] += particle[index].velocity[1] * 0.6;

    /* Applying propellsion to player particle */
    particle[AI].mass -= particle[index].mass;
    particle[AI].radius = sqrt(particle[AI].mass);
    particle[AI].velocity[0] -= particle[index].velocity[0] * particle[PLAYER].mass * dampening;
    particle[AI].velocity[1] -= particle[index].velocity[1] * particle[PLAYER].mass * dampening;
}


void propelAITowards(Real objX, Real objY) {
    int index = numParticles;
    float acceleration = 1.2;
    float dampening = 0.5;
    Real vecX, vecY;

    /* Adds ejected mass to array */
    addParticle(false);
    particle[index].mass = 0.0005;
    particle[index].radius = sqrt(particle[index].mass);
    particle[index].position[0] =  particle[AI].position[0];
    particle[index].position[1] =  particle[AI].position[1];

    /* Vector from AI to click */
    vecX = particle[AI].position[0] - objX;
    vecY = particle[AI].position[1] - objY;

    particle[index].velocity[0] = vecX * acceleration;
    particle[index].velocity[1] = vecY * acceleration;

    /* Moves mass along path - no collision with player */
    particle[index].position[0] += particle[index].velocity[0] * 0.6;
    particle[index].position[1] += particle[index].velocity[1] * 0.6;

    /* Applying propellsion to player particle */
    particle[AI].mass -= particle[index].mass;
    particle[AI].radius = sqrt(particle[AI].mass);
    particle[AI].velocity[0] -= particle[index].velocity[0] * particle[PLAYER].mass * dampening;
    particle[AI].velocity[1] -= particle[index].velocity[1] * particle[PLAYER].mass * dampening;
}


void moveAI() {
    Real sum_radii, sum_radii_sq, n[2], n_mag_sq;
    /* "Bounding" sphere */
    float boundRadius = particle[AI].radius + BOUND;
    bool coll = false;
    
    /* Check against the player particle */
    sum_radii = particle[PLAYER].radius + boundRadius;
    sum_radii_sq = sum_radii * sum_radii;
    n[0] = particle[AI].position[0] - particle[PLAYER].position[0];
    n[1] = particle[AI].position[1] - particle[PLAYER].position[1];
    n_mag_sq = n[0] * n[0] + n[1] * n[1];
    
    if(n_mag_sq <= sum_radii_sq) {
        coll = true;
    }
    
    /* Mote is within range of player, decide to attack or run */
    if(coll) {
        if(particle[AI].radius <= particle[PLAYER].radius) {
            propelAIAway(particle[PLAYER].position[0], particle[PLAYER].position[1]);
            cooling = true;
        }
        else {
            propelAITowards(particle[PLAYER].position[0], particle[PLAYER].position[1]);
            cooling = true;
        }
    }
}


/*******************************************************************
*
*   Collision functions
*
*******************************************************************/

float sumKineticEnergy() {
    Real v_sq, K;

    K = 0;
    for (int i = 0; i < numParticles; i++) {
        v_sq = particle[i].velocity[0] * particle[i].velocity[0] +
                    particle[i].velocity[1] * particle[i].velocity[1];
        K += 0.5 * particle[i].mass * v_sq;
    }

    return K;
}


void sumMomentum(Real *p) {
    p[0] = p[1] = 0;
    
    for (int i = 0; i < numParticles; i++) {
        p[0] += particle[i].mass * particle[i].velocity[0];
        p[1] += particle[i].mass * particle[i].velocity[1];
    }
    
    p[0] += arena.momentum[0];
    p[1] += arena.momentum[1];
}


void eulerStepSingleParticle(Particle &p, Real h) {
    p.position[0] += h * p.velocity[0];
    p.position[1] += h * p.velocity[1];
}


void integrateMotionParticles(Real h) {
    for (int i = 0; i < numParticles; i++) {
        eulerStepSingleParticle(particle[i], h);
    }
}


void collideParticleWall(Particle &p, Arena &a, int index) {
    double pi = 4.0 * atan(1.0);
    float radius = 0;
    Real angle, length, dotNV, outOfRange;
    float dp[2] = {0.0,0.0};
    float prevPos[2] = {0.0,0.0};
    Real normal[2] = {0.0,0.0};
    int count = 0;

    /* Radius at particle's current position */
    radius = pow(p.position[0], 2) + pow(p.position[1], 2);
    radius = sqrt(radius);

    if(radius < arena.max[0] - p.radius - 0.001) {
        return;
    }

    /* Collision */
    prevPos[0] = p.position[0];
    prevPos[1] = p.position[1];

    /* Keeping particle within arena */
    outOfRange = radius + p.radius - a.max[0];
    while(outOfRange >= 0) {
        /* Step back */
        p.position[0] -= arena.dt / 3 * p.velocity[0];
        p.position[1] -= arena.dt / 3 * p.velocity[1];

        radius = pow(p.position[0], 2) + pow(p.position[1], 2);
        radius = sqrt(radius);

        outOfRange = radius + p.radius - a.max[0];
    }

    normal[0] = -prevPos[0];
    normal[1] = -prevPos[1];

    length = sqrt(pow(normal[0], 2) + pow(normal[1], 2));
    normal[0] /= length;
    normal[1] /= length;

    dotNV = normal[0] * p.velocity[0] + normal[1] * p.velocity[1];

    /* Rebound */
    p.velocity[0] = p.velocity[0]  - 2 * dotNV * normal[0];
    p.velocity[1] = p.velocity[1]  - 2 * dotNV * normal[1];

    dp[0] += p.mass * -2.0 * p.velocity[0];
    dp[1] += p.mass * -2.0 * p.velocity[1];

    arena.momentum[0] += dp[0];
    arena.momentum[1] += dp[1];
}


void collideParticlesWall(void) {
    for (int i = 0; i < numParticles; i++) {
        if (debug[debug_wall]) {
            printf("%d %f %f\n", 
                 i, particle[i].position[0], particle[i].position[1]);
        }
        
        collideParticleWall(particle[i], arena, i);
    }
}


void recalculateVelocity(Particle &biggest, Particle smallest) {
    float totalMass = biggest.mass + smallest.mass;
    biggest.velocity[0] = biggest.mass * biggest.velocity[0] + smallest.mass * smallest.velocity[0];
    biggest.velocity[0] /= totalMass;
    
    
    biggest.velocity[1] = biggest.mass * biggest.velocity[1] + smallest.mass * smallest.velocity[1];
    biggest.velocity[1] /= totalMass;
}
                            

void moteCollision(int p1, int p2) {
    float massDiff = particle[p1].mass - particle[p2].mass;
    int biggest = 0;
    int smallest = 0;
    
    if(massDiff < 0) {
        biggest = p2;
        smallest = p1;
    }
    else {
        biggest = p1;
        smallest = p2;
    }

    particle[smallest].collidedWith = biggest;
    particle[smallest].collided = true;
    particle[smallest].absorbed = true;
    particle[biggest].collided = true;
}



void playerCollision(Particle &p2, int index) {
    float mass_diff = particle[PLAYER].mass - p2.mass;
    
    /* Determines if the mote is bigger than the player */
    if(mass_diff <= 0) {
        particle[PLAYER].absorbed = true;
        particle[PLAYER].collided = true;
        particle[PLAYER].collidedWith = index;
        particle[index].collided = true;
    }
    else {
        p2.collided = true;
        p2.absorbed = true;
        particle[PLAYER].collided = true;
        p2.collidedWith = PLAYER;
    }
}


inline bool collisionDetectionParticles(Particle &p1, Particle &p2) {
    Real sum_radii, sum_radii_sq, n[2], n_mag_sq;

    sum_radii = p1.radius + p2.radius;
    sum_radii_sq = sum_radii * sum_radii;
    n[0] = p2.position[0] - p1.position[0];
    n[1] = p2.position[1] - p1.position[1];
    n_mag_sq = n[0] * n[0] + n[1] * n[1];
    
    if (n_mag_sq <= sum_radii_sq) {
        return true;
    }
    else {
        return false;
    }
}


void collisionReactionParticles1D(Particle &p1, Particle &p2) {
    float m1, m2, v1i, v1f, v2i, v2f;

    m1 = p1.mass;
    m2 = p2.mass;
    v1i = p1.velocity[0];
    v2i = p2.velocity[0];
    v1f = (m1 - m2) / (m1 + m2) * v1i + 2.0 * m2 / (m1 + m2) * v2i;
    v2f = 2.0 * m1 / (m1 + m2) * v1i + (m2 - m1) / (m1 + m2) * v2i;

    p1.velocity[0] = v1f;
    p2.velocity[0] = v2f;
}   

     
void collisionReactionParticles2DprojNormal(Particle &p1, Particle &p2) {
    Real n[2], n_mag;
    Real projnv1, projnv2;
    Real m1, m2, v1i, v2i, v1f, v2f;

    /* Normal vector n between centres. */
    n[0] = p2.position[0] - p1.position[0];
    n[1] = p2.position[1] - p1.position[1];

    /* Normalise n. */
    n_mag = sqrt(n[0] * n[0] + n[1] * n[1]);
    n[0] /= n_mag;
    n[1] /= n_mag;

    /* Vector projection/component/resolute of velocity in n direction. */
    projnv1 = n[0] * p1.velocity[0] + n[1] * p1.velocity[1];
    projnv2 = n[0] * p2.velocity[0] + n[1] * p2.velocity[1];

    /* Use 1D equations to calculate final velocities in n direction. */
    v1i = projnv1;
    v2i = projnv2;
    m1 = p1.mass;
    m2 = p2.mass;
    v1f = (m1 - m2) / (m1 + m2) * v1i + 2.0 * m2 / (m1 + m2) * v2i;
    v2f = 2.0 * m1 / (m1 + m2) * v1i + (m2 - m1) / (m1 + m2) * v2i;

    /* Vector addition to solve for final velocity. */
    p1.velocity[0] = (p1.velocity[0] - v1i * n[0]) + v1f * n[0];
    p1.velocity[1] = (p1.velocity[1] - v1i * n[1]) + v1f * n[1];
    p2.velocity[0] = (p2.velocity[0] - v2i * n[0]) + v2f * n[0];
    p2.velocity[1] = (p2.velocity[1] - v2i * n[1]) + v2f * n[1];
}


void collideParticlesBruteForce(Real h) {
    int i, j;

    for (i = 0; i < numParticles - 1; i++) {
        for (j = i + 1; j < numParticles; j++) {
          if (collisionDetectionParticles(particle[i], particle[j])) {
            if (debug[debug_collideParticlesBruteForce]) {
                printf("collideParticlesBruteForce: collision %d %d\n", i, j);
            }

            /* Take step back. Better approaches possible. */
            eulerStepSingleParticle(particle[i], -arena.dt);       
            eulerStepSingleParticle(particle[j], -arena.dt);       

            if (debug[debug_collideParticlesBruteForce]) {
                printf("velocities before: %f %f %f %f\n", 
                    particle[i].velocity[0], particle[i].velocity[1],
                    particle[j].velocity[0], particle[j].velocity[1]);
            }

            /* Collision */
            if (dimension == 1) {
                collisionReactionParticles1D(particle[i], particle[j]);
            }
            else if (dimension == 2) {
                if(i==PLAYER) {
                    playerCollision(particle[j], j);
                }
                else if (reacCalc == basisChange) {
                    moteCollision(i, j);
                }
                else if (reacCalc == projNormal) {
                    collisionReactionParticles2DprojNormal(particle[i], 
                                                               particle[j]);
                }
                else 
                    panic("collision reaction calculation not specified\n");
            }

            if (debug[debug_collideParticlesBruteForce]) {
                printf("velocities after: %f %f %f %f\n", 
                    particle[i].velocity[0], particle[i].velocity[1],
                    particle[j].velocity[0], particle[j].velocity[1]);
            }

            /* Step forward. */
            eulerStepSingleParticle(particle[i], h);        
            eulerStepSingleParticle(particle[j], h);        

            if (debug[debug_collideParticlesBruteForce]) {
                printf("velocities after: %f %f %f %f\n", 
                    particle[i].velocity[0], particle[i].velocity[1],
                    particle[j].velocity[0], particle[j].velocity[1]);
            }

            /* Check walls. */
            // collideParticleWall(particle[i], arena);
            // collideParticleWall(particle[j], arena);

            if (debug[debug_collideParticlesBruteForce]) {
                printf("velocities after: %f %f %f %f\n", 
                    particle[i].velocity[0], particle[i].velocity[1],
                    particle[j].velocity[0], particle[j].velocity[1]);
            }
         }
       }
    }
}


inline void calcGridIndex(Particle &p, Arena a, Real gridCellSize[2], 
                                      int gridNumCells[2], int index[2])
{
    index[0] = (int)((p.position[0] - a.min[0]) / gridCellSize[0]);
    index[1] = (int)((p.position[1] - a.min[1]) / gridCellSize[1]);

#ifdef DEBUG
    printf("\npos[0]  %0.3f  pos[1] %0.3f\n", p.position[0],p.position[1]);
    printf("a.min[0]  %0.3f  a.min[1] %0.3f\n", a.min[0],a.min[1]);
    printf("gridCellSize[0]  %0.3f  gridCellSize[1] %0.3f\n", gridCellSize[0],gridCellSize[1]);
#endif
    
    if (debug_range_check) {
        if (index[0] < 0 || index[0] > gridNumCells[0] - 1) {
            panic("gridIndex: index out of range 1\n");
        }
        if (index[1] < 0 || index[1] > gridNumCells[1] - 1) {
            panic("gridIndex: index out of range 2\n");
        }   
    }
}


void collideParticlesUniformGrid(Real h) {
    Real gridCellSize[2];
    int **gridCellParticleCount, **gridCellParticleListEnd, *gridCellParticleList;
    int gridNumCells[2], gridSize, gridIndex[2], gridCellParticleListStart;
    int gridIndexMin[2], gridIndexMax[2];
    int i, j, k, s, t, p1, p2, total;

    /* Work out grid dimensions and allocate. */
    gridNumCells[0] = (int)(sqrt(numParticles) + 1);
    gridNumCells[1] = (int)(sqrt(numParticles) + 1);
    gridCellSize[0] = (arena.max[0] - arena.min[0]) / gridNumCells[0];
    gridCellSize[1] = (arena.max[1] - arena.min[1]) / gridNumCells[1];
    gridSize = gridNumCells[0] * gridNumCells[1];

    /* Assumption. */
    for (i = 0; i < numParticles; i++) {
        if (particle[i].radius * 2.0 > gridCellSize[0] ||
          particle[i].radius * 2.0 > gridCellSize[1]) {
            panic("collideParticlesUniformGrid: particle diameter > cellSize\n");
        }
    }
 
    /* Allocate arrays. */
    gridCellParticleCount = (int **)malloc(gridNumCells[0] * sizeof(int *));
    if (gridCellParticleCount == 0) {
        panic("collideParticlesUniformGrid: malloc failed\n");
    }

    gridCellParticleListEnd = (int **)malloc(gridNumCells[0] * sizeof(int *));
    if (gridCellParticleListEnd == 0) {
        panic("collideParticlesUniformGrid: malloc failed\n");
    }
    
    for (i = 0; i < gridNumCells[0]; i++) {
        gridCellParticleCount[i] = (int *)malloc(gridNumCells[1] * sizeof(int));
        if (gridCellParticleCount[i] == 0) {
            panic("collideParticlesUniformGrid: malloc failed\n");
        }
            
        gridCellParticleListEnd[i] = (int *)malloc(gridNumCells[1] * sizeof(int));
        if (gridCellParticleListEnd[i] == 0) {
            panic("collideParticlesUniformGrid: malloc failed\n");
        }
    }
    gridCellParticleList = (int *)malloc(numParticles * sizeof(int));

    /* Initialise grid particle count. */
    for (i = 0; i < gridNumCells[0]; i++) {
        for (j = 0; j < gridNumCells[1]; j++) {
            gridCellParticleCount[i][j] = 0;
        }
    }

    /* Cell counts. */
    for (i = 0; i < numParticles; i++) {
        calcGridIndex(particle[i], arena, gridCellSize, gridNumCells, gridIndex);
        gridCellParticleCount[gridIndex[0]][gridIndex[1]] += 1;
    }

    if (debug[debug_collideParticlesUniformGrid]) {
        printf("collideParticlesUniformGrid: gridCellParticleCount\n");
        for (i = 0; i < gridNumCells[0]; i++) {
            for (j = 0; j < gridNumCells[1]; j++) {
                printf("%d %d %d\n", i, j, gridCellParticleCount[i][j]);
            }
        }
    }

    /* Work out end of cell lists by accumulating cell counts. */
    for (i = 0; i < gridNumCells[0]; i++) {
        for (j = 0; j < gridNumCells[1]; j++) {
            gridCellParticleListEnd[i][j] = 0;
        }
    }
  
    total = 0;
    for (i = 0; i < gridNumCells[0]; i++) {
        for (j = 0; j < gridNumCells[1]; j++) {
            total = total + gridCellParticleCount[i][j];
            gridCellParticleListEnd[i][j] = total - 1;
        }
    }

    if (debug[debug_collideParticlesUniformGrid]) {
        printf("collideParticlesUniformGrid: gridCellParticleListEnd\n");
        for (i = 0; i < gridNumCells[0]; i++) {
            for (j = 0; j < gridNumCells[1]; j++) {
                printf("%d %d %d\n", i, j, gridCellParticleListEnd[i][j]);
            }
        }
    }

    /* Build particle lists. */
    for (i = 0; i < gridNumCells[0]; i++) {
        for (j = 0; j < gridNumCells[1]; j++) {
            gridCellParticleCount[i][j] = 0;
        }
    }

    for (i = 0; i < numParticles; i++) {
        calcGridIndex(particle[i], arena, gridCellSize, gridNumCells, gridIndex);
        
        gridCellParticleList[gridCellParticleListEnd[gridIndex[0]][gridIndex[1]] - 
                                gridCellParticleCount[gridIndex[0]][gridIndex[1]]] = i;
                                
        gridCellParticleCount[gridIndex[0]][gridIndex[1]] += 1;
    }
    
    if (debug[debug_collideParticlesUniformGrid]) {
        printf("collideParticlesUniformGrid: gridCellParticleList\n");
        for (i = 0; i < gridNumCells[0]; i++) {
            for (j = 0; j < gridNumCells[1]; j++) {
                gridCellParticleListStart = 
                gridCellParticleListEnd[i][j] - gridCellParticleCount[i][j] + 1;
                printf("particle list %d %d\n", i, j);
                for (k = gridCellParticleListStart; k < gridCellParticleListEnd[i][j]; k++)
                    printf("%d\n", gridCellParticleList[k]);
                    printf("\n");
            }
        }
    }

    /* Collision detection. */
    for (i = 0; i < numParticles; i++) {
        calcGridIndex(particle[i], arena, gridCellSize, gridNumCells, gridIndex);

        /* Grid index bounds for this particle. */
        gridIndexMin[0] = gridIndex[0] - 1;
        if (gridIndexMin[0] < 0) {
            gridIndexMin[0] = 0;
        }
        
        gridIndexMin[1] = gridIndex[1] - 1;
        if (gridIndexMin[1] < 0) {
            gridIndexMin[1] = 0;
        }
        
        gridIndexMax[0] = gridIndex[0] + 1;
        if (gridIndexMax[0] > gridNumCells[0] - 1) {
            gridIndexMax[0] = gridNumCells[0] - 1;
        }
        
        gridIndexMax[1] = gridIndex[1] + 1;
        if (gridIndexMax[1] > gridNumCells[1] - 1) {
            gridIndexMax[1] = gridNumCells[1] - 1;
        }

        p1 = i;

        for (s = gridIndexMin[0]; s <= gridIndexMax[0]; s++) {
            for (t = gridIndexMin[1]; t <= gridIndexMax[1]; t++) {
                gridCellParticleListStart = gridCellParticleListEnd[s][t] 
                                            - gridCellParticleCount[s][t] + 1;
                                            
                for (j = gridCellParticleListStart; j <= gridCellParticleListEnd[s][t]; j++) {
                    p2 = gridCellParticleList[j];
                    /* Don't test particle against itself. */
                    if (p2 == p1)
                        continue;
                    /* Only test pairs once. */
                    if (p2 < p1)
                        continue;

                    if (debug[debug_collideParticlesUniformGrid])
                        printf("collideParticlesUniformGrid: testing %d %d\n", p1, p2);

                    if (collisionDetectionParticles(particle[p1], particle[p2])) {
                        if (debug[debug_collideParticlesUniformGrid])
                            printf("collision: %d %d\n", p1, p2);

                        /* Take step back. Better approaches possible. */
                        eulerStepSingleParticle(particle[p1], -h);  
                        eulerStepSingleParticle(particle[p2], -h);  

                        if (debug[debug_collideParticlesUniformGrid]) {
                            printf("velocities before: %f %f %f %f\n", 
                            particle[p1].velocity[0], particle[p1].velocity[1],
                            particle[p2].velocity[0], particle[p2].velocity[1]);
                        }

                        /* Collision */
                        if(p2 >= numParticles) {
                            continue;
                        }
                        
                        if (dimension == 1) {
                            collisionReactionParticles1D(particle[p1], particle[p2]);
                        }
                        else if (dimension == 2) {
                            if(i == PLAYER) {
                                playerCollision(particle[p2], p2);
                                gridIndexMax[1]--;
                            }
                            else if (reacCalc == basisChange) {
                                moteCollision(p1, p2);
                            }
                            else if (reacCalc == projNormal) {
                                collisionReactionParticles2DprojNormal(particle[p1], particle[p2]);
                            }
                            else {
                                panic("collision reaction calculation not specified\n");
                            }
                        }

                        if (debug[debug_collideParticlesUniformGrid]) {
                            printf("velocities after: %f %f %f %f\n", 
                            particle[p1].velocity[0], particle[p1].velocity[1],
                            particle[p1].velocity[0], particle[p2].velocity[1]);
                        }

                        /* Step forward. */
                        eulerStepSingleParticle(particle[p1], h);   
                        eulerStepSingleParticle(particle[p2], h);   

                        /* Check walls. */
                        // collideParticleWall(particle[p1], arena);
                        // collideParticleWall(particle[p2], arena);
                    }
                }
            }
        }
    }
 
    /* Free arrays. */
    for (i = 0; i < gridNumCells[0]; i++) {
        free(gridCellParticleCount[i]);
        free(gridCellParticleListEnd[i]);
    }
    
    free(gridCellParticleCount);
    free(gridCellParticleListEnd);
    free(gridCellParticleList);
}


/*******************************************************************
*
*   Display functions
*
*******************************************************************/

void writeText(char *buffer, float x, float y) {
    int len, i;
    
    /* Write text */
    glColor3f(1,1,1);
    glRasterPos2f(x,y);
    len = (int)strlen(buffer);
    for (i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, buffer[i]);
    }
}

void setRenderMode(renderMode rm) {
    if (rm == wire) {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_NORMALIZE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    else if (rm == solid) {
        glShadeModel(GL_SMOOTH);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
}


void changeRenderMode(void) {
    if (renMode == wire) {
        renMode = solid;
    } else {
        renMode = wire;
    }
    
    setRenderMode(renMode);
}


void displayMenu() {
    char buffer[80];

    glBegin(GL_QUADS);
        glColor3f(0.55, 0.22, 0.55);
        glVertex2f(5, 5);
        glVertex2f(5, -5);
        glVertex2f(-5, -5);
        glVertex2f(-5, 5);
    glEnd();

    sprintf(buffer, "Level %d", arena.level);
    writeText(buffer, -1, 1.5);

    switch(arena.outcome) {
        case PLAY:
            sprintf(buffer, "[S] to start.");
            writeText(buffer, -1.5, 0);
            break;

        case WIN:
            if(arena.level == 1) {
                sprintf(buffer, "You WON!! =0)");
                writeText(buffer, -2, 0);
                sprintf(buffer, "[S] to move to next level.");
                writeText(buffer, -3, -1.5);
            } else {
                sprintf(buffer, "You are the conqueror!! =0)");
                writeText(buffer, -3.5, 0);
                sprintf(buffer, "[S] to restart. [Q] to quit.");
                writeText(buffer, -3, -1.5);
            }
            break;

        case PAUSE:
            sprintf(buffer, "[S] to restart.");
            writeText(buffer, -1.5, 0);
            break;

        case LOSE_ABSORBED:
            sprintf(buffer, "You lost!!");
            writeText(buffer, -1.5, 0);
            sprintf(buffer, "[S] to restart. [Q] to quit.");
            writeText(buffer, -3, -1.5);
            break;

        case LOSE_SMALLEST:
            sprintf(buffer, "You lost!!");
            writeText(buffer, -3, 0);
            sprintf(buffer, "Your death is inevitable.");
            writeText(buffer, -3, -1);
            sprintf(buffer, "[S] to restart. [Q] to quit.");
            writeText(buffer, -3, -2.5);
            break;

    }
}


void displayArena(void) {
    float x, y;
    int index  = 0;

    glBegin(GL_LINE_LOOP);
        for(index=0; index<360; index++) {
            x = arena.outline[index][0];
            y = arena.outline[index][1];
            glVertex3f(x, y, 0.0);
        }
    glEnd();
}

  
void displayParticle(Particle *p, float sx, float sy, float sz) {
    glPushMatrix();
        glScalef(sx, sy, sz);
        gluDisk(p->quadric, 0.0, p->radius, p->slices, p->loops);
    glPopMatrix();
}


void displayParticles(void) {
    int i;

    /* Display particles. */
    for (i = 0; i < numParticles; i++) {
        if (debug[debug_particle]) {
          printf ("displayParticles: x %f y %f\n", 
                particle[i].position[0], particle[i].position[1]);
        }
      
        glPushMatrix();
            if(i == PLAYER) {
                glColor3f(0.4,1,0.4);
            }
            else if(particle[i]. mass > particle[PLAYER].mass) {
                if(arena.level == 2 && i == AI) {
                    glColor3f(1,0,0.3);
                } else {
                    glColor3f(1,0,0);
                }
            }
            else {
                if(arena.level == 2 && i == AI) {
                    glColor3f(0,0,1);
                } else {
                    glColor3f(0.21,0.6,1);
                }
            }
            
            glTranslatef(particle[i].position[0], particle[i].position[1], 0.0); 
            displayParticle(&particle[i], 1.0, 1.0, 1.0);
        glPopMatrix();
    }
}


void display(void) {
    char buffer[80];
    GLenum err;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glPushMatrix();
        /* Display particle and arena. */
        glPushMatrix();
            /* Adjust for camera */
            glScalef(camera.scale,camera.scale,0);
            glTranslatef(-camera.pos[0], -camera.pos[1], -camera.pos[2]);   
            
            glColor3f (0.8, 0.8, 0.8);         
            displayArena();
            
            displayParticles();
        glPopMatrix();
        
        /* Display menu */
        switch(arena.outcome) {
            case WIN:
            case LOSE_SMALLEST:
            case LOSE_ABSORBED:
                displayMenu();
                break;

            default:
                if(!go) {
                    displayMenu();
                }
                break;
        }
    glPopMatrix();

    glutSwapBuffers();
    
    /* Check for errors. */
    while ((err = glGetError()) != GL_NO_ERROR)
        printf("%s\n",gluErrorString(err));
}


/*******************************************************************
*
*   Update functions
*
*******************************************************************/

void updateParticles(void) {
    static Real time = 0.0;
    float massDecrement = 0.001;
    int biggest;
    static char buffer[80];

    if(!go) {
        return;
    }

    /* Calculate time increment. */
    arena.dt = elapsedTime - time;
    time = elapsedTime;

    if (debug[debug_time]) {
        printf("updateParticles: time %f %f\n", time, arena.dt);
    }

    /* Compute new positions of particles. */
    integrateMotionParticles(arena.dt);

    /* Collisions against walls. */
    collideParticlesWall();

    for(int i=0; i<numParticles; i++) {
        if(particle[i].absorbed) {
            biggest = particle[i].collidedWith;

            particle[i].mass -= massDecrement;
            particle[i].radius = sqrt(particle[i].mass);

            particle[biggest].mass += massDecrement;
            particle[biggest].radius = sqrt(particle[biggest].mass);
        }

        if(particle[i].mass <= 0) {
            if(i == PLAYER) {
                arena.outcome = LOSE_ABSORBED;
                go = false;
                return;
            }

            biggest = particle[i].collidedWith;
            particle[biggest].collided = false;
            recalculateVelocity(particle[biggest], particle[i]);
            removeParticle(i);
        }

        /* Only checks if the player hasn't collided with a larger mote */
        if(!particle[PLAYER].collided) {
            switch(calculateMoteSize()) {
                case WIN:
                    arena.outcome = WIN;
                    go = false;
                    break;

                case LOSE_SMALLEST:
                    arena.outcome = LOSE_SMALLEST;
                    go = false;
                    break;
            }
        }
    }

    /* Collisions amongst particles. */
    if (CDmethod == bruteForce) {
        collideParticlesBruteForce(arena.dt);
    }
    else if (CDmethod == uniformGrid) {
        collideParticlesUniformGrid(arena.dt);
    }
    else {
        panic("updateParticles: unknown collision detection method\n");
    }

    if (debug[debug_sum_kinetic_energy]) {
        printf("K = %f\n", sumKineticEnergy());
    }
    if (debug[debug_sum_momentum]) {
        Real p[2];
        sumMomentum(p);
        printf("p = %f %f\n", p[0], p[1]);
    }
}


void updateCamera() {
    if(camera.zoom == camera.maxZoom) {
        /* Centres on middle of arena */
        camera.pos[0] = 0.0;
        camera.pos[1] = 0.0;
    }
    else {
        /* Centres camera on player */
        camera.pos[0] = particle[PLAYER].position[0];
        camera.pos[1] = particle[PLAYER].position[1];
    }

    /* Zoom the camera in/out (speed is proportional to current zoom) */
    if (camera.zooming && camera.zoom >= MIN_ZOOM && camera.zoom <= camera.maxZoom) {
        camera.zoom -= camera.zoom * camera.direction * camera.sensitivity * 0.05;

        /* Locks camera's zoom to min or max */
        if(camera.zoom > camera.maxZoom) {
            camera.zoom = camera.maxZoom;
        }
        
        if(camera.zoom < MIN_ZOOM) {
            camera.zoom = MIN_ZOOM;
        }

        /* Adjust scale factor */
        camera.scale = 1 / camera.zoom;
    }
}


void update(void) {
    updateCamera();
    
    /* AI mote for level 2 only */
    if(arena.level == 2) {
        /* If cooldown is not in effect, check to see if the mote can move */
        if(!cooling && go) {
            moveAI();
        }
        if(cooling) {
            coolPeriod++;
            if(coolPeriod == COOLTIME) {
                cooling = false;
                coolPeriod = 0;
            }
        }
    }
    
    if (!go) {
        glutPostRedisplay();
        return;
    }

    elapsedTime = glutGet(GLUT_ELAPSED_TIME) / (Real)milli - startTime;

    updateParticles();

    glutPostRedisplay();
}


void resetGame() {
#ifdef DEBUG
    printf("\n\nreset\n");
#endif

    numParticles = 0;
    elapsedTime = 0;
    startTime = 0;
    initialiseArena();
    initialiseParticlesRandomly();

    if(arena.outcome == WIN) {
        arena.level++;

        /* End of game. */
        if(arena.level > 2) {
            arena.level = 1;
        }
    } 
  
    arena.outcome = PLAY;

#ifdef DEBUG
    for(int i=0; i<numParticles; i++){
        printf("%d pos[0]  %0.3f  pos[1] %0.3f\n",i, particle[i].position[0],particle[i].position[1]);
    }

    printf("level = %d ", arena.level);
    printf("outcome = %d \n", arena.outcome);
    printf("min[0]  %0.3f  max[0] %0.3f\n",arena.min[0],arena.max[0]);
    printf("min[1]  %0.3f  max[1] %0.3f\n",arena.min[1],arena.max[1]);
    printf("momentum[0]  %0.3f  momentum[1] %0.3f\n",arena.momentum[0],arena.momentum[1]);
#endif
}



/*******************************************************************
*
*   Key and Mouse functions
*
*******************************************************************/

void keyboardCB(unsigned char key, int x, int y) {
    switch (key) {
        case 'q':
        case 'Q':
        case 27:
           quit();
           break;
           
        case 'w':
           changeRenderMode();
           break;
           
        case 'd':
          if (CDmethod == uniformGrid) {
              CDmethod = bruteForce;
              printf("bruteforce \n");
          }
          else if (CDmethod == bruteForce) {
              CDmethod = uniformGrid;
              printf("uniform grid \n");
          }
          break;
           
        case 'r':
          if (reacCalc == projNormal)
            reacCalc = basisChange;
          else if (reacCalc == basisChange)
            reacCalc = projNormal;
          break;
           
        case 's':
        case 'S':
          if (!go) {
            if(arena.outcome != PAUSE && arena.outcome != PLAY) {
                resetGame();
                return;
            }

            go = true;
            arena.outcome = PLAY;

            if(elapsedTime == 0) {
                startTime = glutGet(GLUT_ELAPSED_TIME) / (Real)milli;
            }
          }
          break;

        case 'p':
        case 'P':
          if (go) {
            go = false;
            arena.outcome = PAUSE;
          }
          break;
    }
    glutPostRedisplay();
}


void updateKeys(int key, bool state) {  
    switch(key) {
        case GLUT_KEY_UP:
            camera.zooming = state;
            if(camera.direction < 0) {
                camera.direction *= -1;
            }
            break;
            
        case GLUT_KEY_DOWN:
            camera.zooming = state;
            if(camera.direction > 0) {
                camera.direction *= -1;
            }
            break;
            
        default:
            break;
    }
}


void arrowDown(int key, int x, int y) {
    updateKeys(key, true);
}


void arrowUp(int key, int x, int y) {
    updateKeys(key, false);
}


void mouseDown(int button, int state, int x, int y) {
    float objX, objY, objZ = 0;

    if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN || !go) {
        return;
    }

    /* Convert to world coordinates */
    float axis = window.axisMax * camera.zoom;

    objX = (axis * 2) * x / ((float)window.width) - axis;
    objY = (axis * 2) * y / ((float)window.height) - axis;
    objY *= -1;

    objX += camera.pos[0];
    objY += camera.pos[1];
    
    propelPlayer(objX, objY);
}

/*******************************************************************
*
*   Init, Reshaping & Main loop
*
*******************************************************************/

void reshape(int w, int h) {
    float aspect = w / (float)h;

    window.width = w;
    window.height = h;
    window.axisMax = 10.0;
    
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-window.axisMax, window.axisMax, -window.axisMax, 
                    window.axisMax, -window.axisMax, window.axisMax);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


/* Initialises given camera with some default values */
void initCamera(void) {
    camera.fov = 75.0f;
    camera.clipNear = 0.1f;
    camera.clipFar = 200.0f;
    camera.zooming = false;
    camera.direction = 1.0f;
    camera.sensitivity = 0.02f;
    camera.pos[0] = 0.0f;
    camera.pos[1] = 0.0f;
    camera.pos[2] = 5.0f;
    camera.zoom = MIN_ZOOM;
    camera.scale = 1/ MIN_ZOOM;
    camera.maxZoom = (arena.max[0] + 2) / window.axisMax;

#ifdef DEBUG
    printf("maxZoom = %0.4f \n", camera.maxZoom);
    printf("window = %0.4f \n", window.axisMax);
#endif
}


void myInit (void) {
#ifdef DEBUG
    printf("my init \n");
#endif
    
    srand (glutGet(GLUT_ELAPSED_TIME));

    setRenderMode(renMode);

    reshape(600, 600);

    initialiseArena();
    arena.level = 1;
    arena.outcome = PLAY;

    initCamera();
    
    initialiseParticlesRandomly();
}


/* 
 *  Main Loop - Open window with initial window size, title bar, 
 *              RGBA display mode, and handle input events.
 */
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Real-Time Rendering : Assignment 3");

    glutDisplayFunc(display);
    glutIdleFunc(update);
    glutReshapeFunc(reshape);

    glutKeyboardFunc(keyboardCB);
    glutSpecialFunc(arrowDown);
    glutMouseFunc(mouseDown);
    glutSpecialUpFunc(arrowUp);

    myInit();

#ifdef DEBUG
    printf("entering main loop...\n");
#endif

    glutMainLoop();
}
