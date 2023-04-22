#pragma once
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"

#include "Particle.hpp"
#include "Spring.hpp"
#include "BendingSpring.hpp"
#include "RobustCCD.hpp"
#include "QuadTree.hpp"

typedef glm::vec<2, double> vec2;

using namespace std;

/**
 * Implementation of a simple particle system
 * @author kry & Mercier-Aubin
 */
class ParticleSystem  {
    
public:
    string name = "";
    std::vector<Particle*> particles;
    std::vector<Spring*> springs;
    std::vector<BendingSpring*> bendingSprings;
    RobustCCD ccd;

    /** and extra spring we'll use to pull on selected particles */
    bool useMouseSpring = false;
    bool useGravityField = false;
    Particle mouseParticle;
    Particle GravityField; //Particle acting as gravity field
    Spring mouseSpring = Spring(&mouseParticle,&mouseParticle);

    /**
     * Creates an empty particle system
     */
    ParticleSystem(){

        ccd = RobustCCD();
        //mouseSpring.k = 10;
        mouseSpring.k = 1e5;   //Making mouse Spring stupidly strong as "moving around masses"
        GravityField.color = glm::vec3(255, 255, 255);
        mouseSpring.l0 = 0;
    }

    ~ParticleSystem() {
        for (int i = 1; i < springs.size(); ++i) {
            delete springs[i];
        }
        for (int i = 1; i < particles.size(); ++i) {
            delete particles[i];
        }
    }
    
    /**
     * Resets the positions of all particles
     */
    void resetParticles() {
        for ( Particle* p : particles ) {
            p->reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles
     */
    void clearParticles() {
        for (Particle* p : particles) { delete p; }
        particles.clear();
        for (Spring* s : springs) { delete s; }
        springs.clear();
        for (BendingSpring* s : bendingSprings) { delete s; }
        bendingSprings.clear();
    }
    
    /** Elapsed simulation time */
    double time = 0;
        
    /**
     *  Compute nodal forces using gravity and spring contributions
     * This function is modified to compute the newly implemented Barnes-Hut gravitation forces simulating outer space bodies and their gravitational pull
     * on each other
     */
    void computeForce() {

        // reset forces and set gravity based on Barnes-Hut or standard atmospheric

        //If not BarnesHut computation
        if (useGravity) { 
            for (Particle* p : particles) {

                vec2 viscousDampingForce = -viscousDamping * p->mass * p->v;
                p->f = vec2(0, 1) * p->mass * gravity + viscousDampingForce; //This is default gravity when simulating in atmosphere
            }
        }

        //If Barnes-Hut (Smart Gravity)
        else if (BarnesHut) {

            bool Attraction = G > 0;
            Smart_Gravity(0.5, Attraction); //BARNES HUT ENTRY CALL GRAVITY COMPUTATIONS WITH THETA AT 0.5

        }


        //Implement "Slow" Gravity Computation Brute Force O(n^2)
        else if (Naive_Gravity) {
            NaiveGravity();
        }

        //applies spring forces
        for (Spring* s : springs) {
            s->apply();
        }

        if (useMouseSpring) {
            mouseSpring.apply();
        }

        if (useGravityField) {
            ApplyGravityField(50);
        }

        for (BendingSpring* s : bendingSprings) {
            s->apply();
        }
    }


    //Naive O(n^2) gravity computation
    void NaiveGravity() {

        //using for any vec to print
        //double min_dist = INFINITY; 

        for (Particle* p : particles) {

            glm::vec2 viscuousForce = -viscousDamping * p->mass * p->v;

            for (Particle* q : particles) {

                if (p == q) continue; //Ignore particle whose forces are being computed

                //Vector pointing from current mass to mass whose gravitational pull is calculated  *p --> *q, --> = r
                glm::vec2 r = q->p - p->p;

                double dist = glm::length(r); //Distance between p and q



                //THIS REDUCES COMPLEXITY IN A NAIVE WAY
                 
                if (dist > 100) continue; //Skip computations for bodies that are too far away

                
                //Add contribution of force
                if (dist < Threshold) {


                    if (q->pinned) {                //Collision with massive body, crash land and compute next particle
                        //q->mass += p->mass;
                        //remove_Particle(p);
                        break;
                    }

                }

                //Scalar in gravitational computations smoothed out with Plummer Smoothing
                double scalar = (G * p->mass * q->mass) / (pow(dist * dist + epsilon * epsilon, 1.5));

                //double scalar = (G * p->mass * q->mass) / ( (glm::length(r) + epsilon) * (glm::length(r) + epsilon));

                p->f.x += scalar * r.x / dist;
                p->f.y += scalar * r.y / dist;

                p->f += viscuousForce;  //Damping as gravity accelerates particles


                }
        }
    }

    //For two bodies (x1,y1,m1) , (x2,y2,m2)
    //m = m1+m2                             //Total Mass
    //cx = (x1*m1 + x2*m2)/m                //Center of masses
    //cy = (y1*m1 + y2*m2)/m            
   
    
    void Smart_Gravity(double theta, bool Attraction) {

        //Initialize first QuadTree (root) as the whole display

        //Instantiate QuadTree root Node with sizes matching whole display
        QuadTree root = QuadTree(width / 2, height / 2, width, height, Attraction, restitutionValue); //Build QuadTree for whole display


        //For each particles in the system
        for (Particle* particle : particles) {
            root.insert(particle); //Insert the particle in the tree to build it
        }
        

        //From this point, the QuadTree is built

        //For each particle in the system, compute the force acting on it
        for (Particle* particle : particles) {
            root.computeForce(particle, theta, 0.45 ); //epsilon is 0.45 here, theta is defined by Smart_Gravity caller
        }

        //root.collision_resolve(restitutionValue); //CCD step of root ZELDA

    }


    //Boolean function checking if a particle is close enough to be at risk of being
    //out of bounds from the simulation
    bool isOutOfBounds(Particle* p, double safety){
        return p->p.x < safety || p->p.x > width - safety || p->p.y < safety || p->p.y > height - safety;
    }

    void damp(Particle* p, double safety, double damping_factor) {
        // Check if particle is close enough to be at risk of going out of bounds
        if (!isOutOfBounds(p, safety)) return; // If not, return

        // Calculate distance from particle to closest boundary
        double dist_left = p->p.x - safety;
        double dist_right = width - safety - p->p.x;
        double dist_top = p->p.y - safety;
        double dist_bottom = height - safety - p->p.y;

        // Damp velocity and forces based on distance from boundaries
        if (dist_left < 0) {
            p->v.x *= 1 - damping_factor * dist_left;
            p->f.x -= p->v.x * damping_factor * dist_left;

        }
        if (dist_right < 0) {
            p->v.x *= 1 - damping_factor * dist_right;
            p->f.x -= p->v.x * damping_factor * dist_right;
        }
        if (dist_top < 0) {
            p->v.y *= 1 - damping_factor * dist_top;
            p->f.y -= p->v.y * damping_factor * dist_top;
        }
        if (dist_bottom < 0) {
            p->v.y *= 1 - damping_factor * dist_bottom;
            p->f.y -= p->v.y * damping_factor * dist_bottom;
        }

        // Stop particle from leaving the boundary by setting its position and velocity to 0 in that direction
        if (p->p.x <= 0) {
            p->p.x = 0;
            if (p->v.x < 0) {
                p->v.x = 0;
            }
            if (p->f.x < 0) {
                p->f.x = 0;
            }
        }

        if (p->p.x >= width) {
            p->p.x = width;
            if (p->v.x > 0) {
                p->v.x = 0;
            }
            if (p->f.x > 0) {
                p->f.x = 0;
            }
        }
        if (p->p.y <= 0) {
            p->p.y = 0;
            if (p->v.y < 0) {
                p->v.y = 0;
            }
            if (p->f.y < 0) {
                p->f.y = 0;
            }
        }
        if (p->p.y >= height) {
            p->p.y = height;
            if (p->v.y > 0) {
                p->v.y = 0;
            }
            if (p->f.x > 0) {
                p->f.x = 0;
            }
        }
    }

   
    void stepVelocities(float h) {
        for (Particle* p : particles) {
            if (p->pinned) {
                p->f = vec2(0, 0); // just to make sure!
                p->v = vec2(0, 0);
            }
            else {
                p->v += (h / p->mass) * p->f;

                damp(p, 1, 0.1);

            }
        }
    }

  
    /**
     * Update the positions with current velocities given the time step
     * @param h time step
     */
    void stepPositions(double h) {
        for (Particle* p : particles) {
            if (p->pinned) continue;
            p->p += h * p->v;
            p->f = vec2(0, 0);
        }
    }

    /** Time in seconds that was necessary to advance the system */
    float computeTime = 0;
    
    /**
     * Advances the state of the system
     * @param elapsed
     */
    void advanceTime( double elapsed ) {
        for (Spring* s : springs) {
            s->k = springStiffness;
            s->c = springDamping;
        }
        
        double now = glfwGetTime();
        
        // Simple symplectic Euler integration
        computeForce();
        stepVelocities(elapsed);
        ccd.applyRepulsion(elapsed, &particles, &springs); //Removed repulsions for Barnes Hut simulation. Replaced with BVH in QuadTree ZELDA
        //bool resultOK = ccd.checkCollision(elapsed, &particles, &springs);
        stepPositions(elapsed);
        
        time += elapsed;
        computeTime = (glfwGetTime() - now);
    }
    
    void filter(VectorXf& v) {
        for ( Particle* p : particles ) {
            if ( !p->pinned ) continue;
            v[ p->index*2+0] = 0;
            v[ p->index*2+1] = 0;
        }
    }

    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    Particle* createParticle( float x, float y, float vx, float vy ) {
        Particle* p = new Particle( x, y, vx, vy );
        p->index = particles.size();
        particles.push_back( p );
        return p;
    }
    
    void remove( Particle* p ) {
    	for ( Spring* s : p->springs ) {
    		Particle* other = s->p1 == p ? s->p2 : s->p1; 
            // remove from the other particle's spring list
            vector<Spring*>::iterator osi = std::remove(other->springs.begin(), other->springs.end(), s);
            other->springs.erase(osi);
            // and remove from the master list
            vector<Spring*>::iterator si = std::remove(springs.begin(), springs.end(), s);
            springs.erase(si);
            delete s;
    	}

    	particles.erase( std::remove(particles.begin(),particles.end(), p ) );
    	// reset indices of each particle :(
    	for ( int i = 0 ; i < particles.size(); i++ ) {
    		particles[i]->index = i;
    	}
    }

    //Test to see if it correctly removes particle
    void remove_Particle(Particle* p) {
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            if (*it == p) {
                particles.erase(it);
                break;
            }
        }
    }

    void toggle_G(bool Attraction) {

        if (Attraction) G = abs(G); //Positive G
        else G = -abs(G); //Negative G
    }

    //Getter method to get constant being used 
    double get_G() {
        return G;
    }

    void ApplyGravityField(double distance) {

        for (Particle* p : particles) {

            //If close enough to Field
            double r = GravityField.distance(p->p.x, p->p.y);

            glm::vec2 r_vec = p->p - GravityField.p;

            double force = (G < 0) ? -1 : 1;

         


            if (r < distance && r > 5) {

                //Apply Force of Field onto p

                force *= 10 * p->mass * p->mass / (r * r + epsilon);  

                //double force = 3 * particle->mass * mass / pow(distance * distance + eps * eps, 1.5);

                //Approximate gravitational forces acting on particle by the center of mass of the QuadTree if its center of mass is far enough (determined by theta)

                p->f.x += force * r_vec.x / r;     //Add the forces to the current particle appropriately
                p->f.y += force * r_vec.y / r;


                if (r < 5.1){            //When close, kill forces and velocities
                    p->f.x *= 0.9;
                    p->f.y *= 0.9;
                    p->v.x *= 0.9;
                    p->v.y *= 0.9;
                }
            }
        }
    }

    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    Spring* createSpring( Particle* p1, Particle* p2 ) {
        Spring* s = new Spring( p1, p2 ); 
        springs.push_back( s );         
        return s;
    }
    
    /**
     * Removes a spring between p1 and p2 if it exists, does nothing otherwise
     * @param p1
     * @param p2
     * @return true if the spring was found and removed
     */
    bool removeSpring( Particle* p1 = new Particle, Particle* p2 = new Particle ) {
    	Spring* found = NULL;
    	for ( Spring* s : springs ) {
    		if ( ( s->p1 == p1 && s->p2 == p2 ) || ( s->p1 == p2 && s->p2 == p1 ) ) {
    			found = s;
    			break;
    		}
    	}
    	if ( found != NULL ) {
            found->p1->springs.erase(std::remove(found->p1->springs.begin(), found->p1->springs.end(), found));
            found->p2->springs.erase(std::remove(found->p2->springs.begin(), found->p2->springs.end(), found));
            springs.erase(std::remove(springs.begin(),springs.end(), found));
			return true;
    	}
    	return false;
    }
    
    void init() {
        // do nothing
    }

    int height = 0;
    int width = 0;

    void display() {

        glPointSize( 4 );
        glBegin( GL_POINTS );
        for ( Particle* p : particles ) {
            double alpha = 0.75;
            if ( p->pinned && p != &GravityField ) {
                glColor4d( 1, 0, 0, alpha );
            } else {
                glColor4d( p->color.x, p->color.y, p->color.z, alpha );
            }
            glVertex2d( p->p.x, p->p.y );
        }

        glEnd();

        glColor4d(0,.5,.5,1);
        glLineWidth(2.0f);
        glBegin( GL_LINES );
        for (Spring* s : springs) {
            glVertex2d( s->p1->p.x, s->p1->p.y );
            glVertex2d( s->p2->p.x, s->p2->p.y );
        }
        if (useMouseSpring) {
            glColor4d(1, 0, 0, 1);
            Spring* s = &mouseSpring;
            glVertex2d(s->p1->p.x, s->p1->p.y);
            glVertex2d(s->p2->p.x, s->p2->p.y);
        }
        glEnd();
    }
    
    bool useGravity = false;
    bool BarnesHut = true; //Added BarnesHut gravity boolean for particle system 
    bool Naive_Gravity = false; //Boolean for Naive Gravity
    double max_dist = 50; //Distance in Naive Gravity where neighbouring particles are ignored
     
    //The lower this is, the less accurate as a whole the system is,but aggregates of particles seem to be following Newton's Gravitational Law pretty well
    //System becomes more dynamic locally
    

    double bound_thres = 4.5; //Distance from boundaries where we start projecting velocities along boundary
    //double G = 6.6743e-11; //For BarnesHut Gravitational Constant
    double G = 3; //Scaled G as masses are already scaled down so their proportional product in Gravity computation looks normal
    double Threshold = 5.5;
    //double Threshold = 4.5; //Distance threshold defining if a body is closer to a massive body by that distance, it will fall into it (crash on it)
    double epsilon = 0.45; //Fraction of distance when particles are "in contact" that serves as correction preventing division by zero
                                      //When computing gravitational forces between particles (This prevents forces from "exploding" upon contact

    double restitutionValue = 0.2; //Restitution coefficient for CCD through QuadTree


    //double const bound_scale = tan(0.9 * M_PI / 2) / Threshold; //This constant is used in a function f(x) = (2/pi) * arctan(kx)
                                                                //Where k is the bound_scale. This function is used to damp particles
                                                                //Approaching the Threshold


    double gravity = 9.8;
    double springStiffness = 100;
    double springDamping = 0;
    double viscousDamping = 0;
    double bendingStiffness = 1e3;
};

