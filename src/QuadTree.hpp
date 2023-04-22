#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "ParticleSystem.hpp"

using namespace std;


//Quadtree Definition For Barnes-Hut Computations
class QuadTree {
private:

    double G;
    static const int capacity = 4;  // maximum number of particles in a quadrant
    bool divided = false;           // flag indicating if the tree is divided
    double x, y, size_x, size_y, size;              // coordinates and size of the quadrant (x,y,diagonal)
    double restitutionValue;
    bool Attraction;
    std::vector<Particle*> particles;  // list of particles in the quadrant (DO I NEED THIS OR CAN I USE THE ONE IN PARTICLE SYSTEM?)

    //std::vector<std::pair<Particle*, Particle*>> particlePairs; //Vector containing pairs of all potential collisions



    QuadTree* nw = nullptr;  // pointer to northwest subquadrant
    QuadTree* ne = nullptr;  // pointer to northeast subquadrant
    QuadTree* sw = nullptr;  // pointer to southwest subquadrant
    QuadTree* se = nullptr;  // pointer to southeast subquadrant

public:

   

    QuadTree(double x, double y, double size_x, double size_y, bool Attraction, double restitutionValue) : x(x), y(y), size_x(size_x), size_y(size_y), Attraction(Attraction), restitutionValue(restitutionValue){
        size = sqrt(size_x * size_x + size_y * size_y); //set size of diagonal based on length of sides
    }


    //Constructor of a QuadTree Quadrant (Node), x and y are center coordinates
    //size_x and size_y are dimensions in x and y (width and height) of the quadrant


    // check if a particle is inside the quadrant
    bool contains(Particle* particle) const {
        return particle->p.x >= x - size_x / 2 && particle->p.x < x + size_x / 2 && //Check if x,y coordinates of particle are bounded by Quadrant coordinates
            particle->p.y >= y - size_y / 2 && particle->p.y < y + size_y / 2;      //As x,y are center, bounds extend left and right by size/2 of corresponding dimension
    }                                                                               //p->p.x between x +- size_x/2, same for y

    // insert a particle into the tree
    bool insert(Particle* particle) {   //If particle does not belong in the quadrant of the tree, return false (couldn't insert it in that QuadTree (Node) )
        if (!contains(particle)) {
            return false;
        }
        if (particles.size() < capacity) {  //If particle does belong in the quadrant and there is capacity left in it, we add said particle to the particles array
            particles.push_back(particle);  //of the quadrant and return true
            return true;
        }
        if (!divided) {                     //If the particle belongs to the quadrant, but there is no more space in said quadrant and the quadrant
            subdivide();                    //has not been divided yet, we subdivide the currentQuadTree into four smaller quadrants hoping these four 
                                            //will be able to each hold less than "capacity" particles.

           
 
        }                                   

        if (nw->insert(particle)) return true; //If the quadrant has already been divided, we recursively try to insert the current particles
        if (ne->insert(particle)) return true; //In the appropriate sub-quadrant of the current quadrant until it reaches its "final" quadrant in the QuadTree.
        if (sw->insert(particle)) return true; //return true if it finds the appropriate quadrant. Return false if the particle is "out of bounds" and couldn't fit in
        if (se->insert(particle)) return true; //the QuadTree (shouldn't happen but preventive).
        return false;
    }


    //nw    ne          //Quadrants are separated like this,  x grows towards the right, y grows downwards, origin (0,0) is located at top-left of the screen
    //sw    se

    // subdivide the quadrant into four subquadrants
    void subdivide() {
        nw = new QuadTree(x - size_x / 4, y - size_y / 4, size_x / 2, size_y / 2, Attraction, restitutionValue);  //Make appropriate constructors with new centers and sizes
        ne = new QuadTree(x + size_x / 4, y - size_y / 4, size_x / 2, size_y / 2, Attraction, restitutionValue);
        sw = new QuadTree(x - size_x / 4, y + size_y / 4, size_x / 2, size_y / 2, Attraction, restitutionValue);
        se = new QuadTree(x + size_x / 4, y + size_y / 4, size_x / 2, size_y / 2, Attraction, restitutionValue);
        divided = true;                                                             //Flag current QuadTree (Node) as divided
    }

    // compute the net gravitational force acting on a particle
    void computeForce(Particle* particle, double theta, double eps) {
        // Compute the distance between the particle and the center of mass of this quadrant

        double centerOfMassX = 0, centerOfMassY = 0, mass = 0;
        for (Particle* p : particles) {                                 //For each particles belonging to the quadrant. (This is done once the QuadTree has been built)
            centerOfMassX += p->p.x * p->mass;                          //Which means that all particles considered are correctly placed in the QuadTree
            centerOfMassY += p->p.y * p->mass;
            mass += p->mass;                                            //Compute the coordinates x,y of its center of mass as well as its "total" mass
        }
        centerOfMassX /= mass;                                          //Center of mass computed
        centerOfMassY /= mass;
        double dx = centerOfMassX - particle->p.x;                      //Compute distance of said "approximated body" to the particle whose forces are being computed
        double dy = centerOfMassY - particle->p.y;
        double distance = sqrt(dx * dx + dy * dy);                      //Get norm of vector (positive value) of distance to approximated body




        // If this quadrant is far enough away, use the center of mass to approximate the force. This uses a ratio of the diagonal of a quadrant to the distance of
        //The current particle to the center of mass of another quadrant
        if (size / distance < theta) {
            //double force = G * particle->getMass() * mass / (distance * distance + eps);  //TO DO get G correctly


            double force = 3 * particle->mass * mass / (distance * distance + eps);  //TO DO get G correctly

            //double force = 3 * particle->mass * mass / pow(distance * distance + eps * eps, 1.5);

            //Approximate gravitational forces acting on particle by the center of mass of the QuadTree if its center of mass is far enough (determined by theta)


            if (!Attraction) force = -force; //Invert force direction if not attraction

            particle->f.x += force * dx / distance;     //Add the forces to the current particle appropriately
            particle->f.y += force * dy / distance;
        }
        // Otherwise, recurse on the child quadrants
        else {
            if (nw) nw->computeForce(particle, theta, eps);
            if (ne) ne->computeForce(particle, theta, eps);
            if (sw) sw->computeForce(particle, theta, eps);
            if (se) se->computeForce(particle, theta, eps);
        }
    }

    bool collision(Particle* p1, Particle* p2, double restitutionValue) {        //Resolve collision between two particles

        if (p1 == p2) return false; //No collision

        glm::vec2 diff = p2->p - p1->p;

        double dist = glm::length(diff); //distance between two particle

        if (dist > 5) return false; //Too far for collision

        glm::vec2 rel_Vel = p2->v - p1->v; //relative velocity between particles

        if (glm::dot(rel_Vel, diff) < 0) { //If not, then particles are moving away from each other

            float impulseMag = -(1 + restitutionValue) * glm::dot(rel_Vel, diff) / (glm::dot(diff, diff) * (1 / p1->mass + 1 / p2->mass));

           // cout << impulseMag;

            p1->v.x -= impulseMag / p1->mass * diff.x;
            p1->v.y -= impulseMag / p1->mass * diff.y;

            p2->v.x += impulseMag / p2->mass * diff.x;
            p2->v.y += impulseMag / p2->mass * diff.y;

            p1->pinned = true;
            p2->pinned = true; //ZELDA

            return true;

        }

        return false;
    }

    //Function that resolves all collision in each leaf node of the quadtree. THIS IS NOT FUNCTIONAL YET
    void collision_resolve(double restitutionValue) {

        bool collided = true;
        int iterations = 0;

     
        if (!divided) { //If leaf - node

            while (collided && iterations++ < 10){
                collided = false;
                for (Particle* p : particles) {
                    for (Particle* q : particles) {
                        if (collision(p, q, restitutionValue)) collided = true; //Collision resolve between two particles in sub-quadrant
                    }
                }
            }
             //Done with current leaf node
        }

        else {  //If not leaf node, iterate over every children until leaf is reached

            if (nw) nw->collision_resolve(restitutionValue);
            if (ne) ne->collision_resolve(restitutionValue);
            if (sw) sw->collision_resolve(restitutionValue);
            if (se) se->collision_resolve(restitutionValue);
        }
    }


    ~QuadTree() {       //Delete subtrees
        delete nw;
        delete ne;
        delete sw;
        delete se;
    }
    
};



//##########################################
//IGNORE THE FOLLOWING LINES, I HAD TROUBLE SETTING THE PATH OF MY QUADTREE.HPP FILE SO I TOOK THE IMAGERECORDER ONE AND RENAMED IT 
//##########################################



//class ImageRecorder {
//public:

	//string path = "../out/";
	//int frameNumber = 0;
	//string filename = "";
	//vector<unsigned char> pixels;
	//vector<unsigned char> rowData;

	//void writeCurrentFrameToFile(GLFWwindow* window) {
		//// get the image
		//int comp = 3; // RGB
		//int width, height;
		//glfwGetFramebufferSize(window, &width, &height);
		//if (pixels.size() != width * height * comp) {
		//	pixels.resize(width * height * comp);
		//	rowData.resize(width * comp);
		//}
		//glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0]);
		//// gross, need to flip the buffer
		//for (int i = 0; i < height / 2; i++) {
			//int r1 = i * width * comp;;
			//int r2 = (height - 1 - i) * width * comp;
			//memcpy(&rowData[0], &pixels[r1], rowData.size());
			//memcpy(&pixels[r1], &pixels[r2], rowData.size());
			//memcpy(&pixels[r2], &rowData[0], rowData.size());
		//}
		//// make the filename
		//std::ostringstream oss;
		//oss << path << "image" << std::setfill('0') << std::setw(5) << frameNumber++ << ".png";
		//filename = oss.str();
		//// write to file
		//int stride_in_bytes = width * comp * sizeof(unsigned char);
		//int rc = stbi_write_png(filename.c_str(), width, height, comp, &pixels[0], stride_in_bytes);
		//if (rc == 0) {
			//cout << "Couldn't write to " << filename << endl;
		//}
	//}
//};

