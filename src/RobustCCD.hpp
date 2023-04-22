#pragma once
#include "Particle.hpp"
#include "Spring.hpp"

class RobustCCD {
public:
	
	bool collision;
	bool useRepulsion;
	int maxIterations;
	double restitutionValue;
	double thresholdDistance;
	double eps;
	double repulsionStiffness;

	/** TODO: number of iterations in the last CCD processing loop, to keep an eye on how tricky the problem is */
	int iterations;
	/** TODO: keep track of highest count needed so far */
	int highestIterations = 0;

	RobustCCD::RobustCCD() :
		eps(0.000001),
		collision(true),
		useRepulsion(true),
		maxIterations(50),
		restitutionValue(0.5),
		thresholdDistance(2),
		iterations(0),
		repulsionStiffness(100) {}

	/**
	 * Try to deal with contacts before they happen
	 * @param h
	 * @param system
	 */
	void applyRepulsion(double h, std::vector<Particle*>* particles, std::vector<Spring*>* springs) {
		//TODO: apply repulsion on all particles

		// - check for point-edge proximity between all particle-edge pairs
		// - find the normal
		// - take care to deal with segment end points carefully
		// - compute an appropriate  impulse if the distance is less than H
		// - make sure your impulse is in the correct direction!
		// - don't apply the impulse if the relative velocity is separating fast enough (see Bridson et al. 2002)
		// - distribute impulse to the three particles involved in the appropriate manner

		if (!useRepulsion) return;


		for (Particle* p1 : *particles) {
			for (Particle* p2 : *particles) {
				glm::vec2 diff = p2->p - p1->p;
				double dist = glm::length(diff); //distance between two particles

				if (p1 == p2 || dist > 4.5) continue; //Not close enough or same particle


				glm::vec2 rel_Vel = p2->v - p1->v; //relative velocity between particles

				if (glm::dot(rel_Vel, diff) < 0) {

					float impulseMag = -(1 + restitutionValue) * glm::dot(rel_Vel, diff) / (glm::dot(diff, diff) * (1 / p1->mass + 1 / p2->mass));
				


					p1->v.x -= impulseMag / p1->mass * diff.x;
					p1->v.y -= impulseMag / p1->mass * diff.y;


					p2->v.x += impulseMag / p2->mass * diff.x;
					p2->v.y += impulseMag / p2->mass * diff.y;

					collision = true;	
				}
			}
		}

		for (Spring* s : *springs) {

			for (Particle* p0 : *particles) {

				if (p0 == s->p1 || p0 == s->p2) {
					continue;

				}


				glm::vec2 diff_AC = s->p1->p - p0->p;
				glm::vec2 diff_AB = s->p1->p - s->p2->p;

				double alpha = (glm::dot(diff_AC, diff_AB)) / (diff_AB.x * diff_AB.x + diff_AB.y * diff_AB.y);

				
				glm::vec2 midp1 = {
					alpha * diff_AB.x,
					alpha * diff_AB.y,

				};

				//p1p0 = diff_AC

				glm::vec2 normal = diff_AC + diff_AB;

				if (alpha <= 0) {
					normal = p0->p - s->p1->p;

				}

				else if (alpha >= 1) {
					normal = p0->p - s->p2->p;

				}

				float d_test = thresholdDistance - glm::length(normal);

				float d = fabs(d_test);


				glm::vec2 Vminus{
					alpha * s->p1->v.x + (1 - alpha) * s->p2->v.x - p0->v.x,
					alpha * s->p1->v.y + (1 - alpha) * s->p2->v.y - p0->v.y,
				};


				float Vn = glm::dot(normal, Vminus);

				//if naturally moving away do nothing
				if (Vn >= 0.1f*d / h) continue;
				//If not withing threshold do nothing
				if (d < -eps || d >= thresholdDistance) continue;


				//Impulse divided by mass as we're only updating velocities which will cancel the term
				double I_c = Vn / 2;


				double vel_check = 0.1 * d / h;

				double p0_collision_mass = p0->mass;
				double p1_collision_mass = s->p1->mass;
				double p2_collision_mass = s->p2->mass;


				if (p0->pinned) {
					p0_collision_mass = INFINITY;
				}

				if (s->p1->pinned) {
					p1_collision_mass = INFINITY;
				}

				if (s->p2->pinned) {
					p2_collision_mass = INFINITY;
				}


				double Ir0 = -min(h * repulsionStiffness * d, p0_collision_mass * (vel_check - Vn));
				double Ir1 = -min(h * repulsionStiffness * d, p1_collision_mass * (vel_check - Vn));
				double Ir2 = -min(h * repulsionStiffness * d, p2_collision_mass * (vel_check - Vn));


				//Apply inelastic and repulsion impulses

				p0->v.x += normal.x * (I_c + 0.25*Ir0);
				p0->v.y += normal.y * (I_c + 0.25*Ir0);

				s->p1->v.x -= alpha * (I_c + 0.25*Ir1) * normal.x;
				s->p1->v.y -= alpha * (I_c + 0.25*Ir1) * normal.y;

				s->p2->v.x -= (1 - alpha) * (I_c + 0.25 * Ir2) * normal.x;
				s->p2->v.y -= (1 - alpha) * (I_c + 0.25 * Ir2) * normal.y;

	
			}
		}

	}

	/**
	* Checks all collisions in interval t to t+h
	* @param h
	* @param system
	* @return true if all collisions resolved
	*/
	bool checkCollision(double h, std::vector<Particle*>* particles, std::vector<Spring*>* springs) {
		// TODO: Most of your assignment code will go here
		// - put your iterative solving code here!
		// - you can use the nextCollision function below to deal with an individual particle-edge pair

		// For each particle-edge pair, find the roots for when the three particles are
		// co-linear, and then pick the first root on (0,h] which corresponds to an actually collision.
		// compute a collision response.  That is, compute an appropriate collision normal, compute the 
		// impulse, and then apply the impulse to the associated particles.  Be sure to deal with pinning 
		// constraints!



		//Initiate first iterative loop iteration for collision detection in interval
		iterations = 0;

		//Set collision flag to true at each time step
		collision = true;

		//Iterative solve
		while (collision && iterations < maxIterations) {

			//Assume no collision in current iteration
			collision = false;
			for (Spring* s : *springs) {

				for (Particle* p : *particles) {

					//For each springs, for each particles

					//Checking collision between spring and all other particles, except for the particles part of the spring
					if (p == s->p1 || p == s->p2) {
						continue;
					}

					//Check if there is a collision

					//At least one collision processed in time step
					if (nextCollision(h, restitutionValue, p, s->p1, s->p2)) collision = true;


				}
			}

			//If no collisions, collision is false and while loop exits, else collision is true (at least one collision processed in (0,h]
			//Iterate again over same interval and update iteration counter
			iterations += 1;
			//cout << iterations;
		}

		if (iterations >= maxIterations) return false;



		return true;
	}

	/**
	* Processes next collision on (0,h] if it exists, between p0 and segment connecting p1 and p2.
	*
	* Watch out for
	* - normal direction
	* - epsilon tests
	* - segment connecting p1 and p2 at zero length! (unlikely?)
	*
	* @param h				timestep
	* @param restitution	bouncyness parameter
	* @param p0				particle
	* @param p1				line segment end point
	* @param p2				line segment end point
	* @return true if collision processed
	*/
	bool RobustCCD::nextCollision(double h, double restitution, Particle* p0, Particle* p1, Particle* p2) {
		//TODO: 
		// * - finds roots
		// * - checks that p0 actually falls on segment
		// * - processes collision by computing and applying impulses
		// * - returns true if collision processed


		//Quadratic equation to solve for collision of the form at^2 + bt + c = 0

		double A = p1->p.x;
		double B = p1->v.x;
		double C = p1->p.y;
		double D = p1->v.y;

		double E = p2->p.x;
		double F = p2->v.x;
		double G = p2->p.y;
		double H = p2->v.y;

		double I = p0->p.x;
		double J = p0->v.x;
		double K = p0->p.y;
		double L = p0->v.y;

		double a = (B * H - B * L - D * F + D * J + F * L - H * J);
		double b = (A * H - A * L + B * G - B * K - C * F + C * J - D * E + D * I + E * L + F * K - G * J - H * I);
		double c = (A * G - A * K - C * E + C * I + E * K - G * I);



		std::vector<double> root;
		bool bound = false;
		double boundary;

		//Linear case
		if (fabs(a) < eps) {

			if (fabs(b) < eps) {

				if (c > -eps && c <= h + eps) root.push_back(0.0);

		
			}

			else {
				//cout << "a = 0";
				root.push_back(-c / b);

			}
		}

		//True quadratic
		else {

			//cout << "quadratic";

			double discriminant = b * b - 4 * a * c;

			if (fabs(discriminant) < eps) {
				root.push_back(-b / (2 * a));
			}

			else if (discriminant >= eps) {
				root.push_back((-b - sqrt(discriminant)) / (2 * a));
				root.push_back((-b + sqrt(discriminant)) / (2 * a));
			}
		}




		//if (fabs(boundary) > eps) continue;

		
		//For each valid root found
		for (double t : root) {

			//if root is valid in time step with epsilon check
			if (t > -eps && t <= h + eps ) {

				
				//Update position of each particles to collision point
				
				glm::vec2 p0_step = p0->p + t * p0->v;
				glm::vec2 p1_step = p1->p + t * p1->v;
				glm::vec2 p2_step = p2->p + t * p2->v;

				//initialize vectors between particles
				
				glm::vec2 diff_AC = p1_step - p0_step;
				glm::vec2 diff_AB = p1_step - p2_step;

			
				// A -> p1
				// B -> p2
				// C -> p0

				
				double alpha = (glm::dot(diff_AC, diff_AB)) / (diff_AB.x * diff_AB.x + diff_AB.y * diff_AB.y);
					
				if (alpha >= -eps && alpha <= 1 + eps) {
					//Root and in segment
					//compute impulses

					//compute normal pointing from segment AB towards particle C

					glm::vec2 AB_normalized = glm::normalize(diff_AB);

					double dot_product = glm::dot(AB_normalized, diff_AC);

					glm::vec2 normal;

					if (dot_product > 0) {
						normal = glm::vec2(-AB_normalized.y, AB_normalized.x);
					}

					else {
						normal = glm::vec2(AB_normalized.y, -AB_normalized.x);

					}

					//relative velocity before collision

					glm::vec2 prev_vel = alpha * p1->v + (1 - alpha) * p2->v - p0->v;

					//If pinned particles, assumption of infinite mass
					double p0_collision_mass = p0->mass;
					double p1_collision_mass = p1->mass;
					double p2_collision_mass = p2->mass;

					if (p0->pinned) {
						p0_collision_mass = INFINITY;
					}

					if (p1->pinned) {
						p1_collision_mass = INFINITY;
					}

					if (p2->pinned) {
						p2_collision_mass = INFINITY;
					}

				
					double j = (1.0 +  restitution) * glm::dot(normal, prev_vel) / ( ((alpha*alpha)/p1_collision_mass) + ((1-alpha) * (1-alpha) / (p2_collision_mass)) + (1/(p0_collision_mass))  );

					//Update p0 velocity

					//p0->v += p0->f * t / (p0->mass);			//Velocity at predicted collision time
					p0->v.x += normal.x * (j  / p0->mass);		//Velocity updated by impulses
					p0->v.y += normal.y * (j  / p0->mass);

					//Update p1 velocity
					//p1->v += p1->f * t / (p1->mass);
					p1->v.x -= alpha*normal.x * (j / p1->mass);		//Velocity updated by impulses
					p1->v.y -= alpha*normal.y * (j / p1->mass);


					//Update p2 velocity
					//p2->v += p2->f * t / (p2->mass);
					p2->v.x -= (1-alpha)*normal.x * (j / p2->mass);		//Velocity updated by impulses
					p2->v.y -= (1 - alpha) * normal.y * (j / p2->mass);



					//p0->pinned = true;
					//p1->pinned = true;
					//p2->pinned = true;


					return true;

				}

			}



		}

		return false;

	}

};

