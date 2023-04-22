Yajiv Luckheenarain 

This is a Barnes-Hut algorithm Implementation to solve N-Body Simulations. This code has been written using a template code by
@paul kry, my COMP559 professor at McGill University.


The simulation is meant to display the difference in performance between a Naive Brute Force Approach at solving the 
gravitational forces between N particles (O(n^2)) and the Barnes-Hut Algorithm O(nlog(n)).

TestSystems are already implemented to experiment with different systems. 

Particles have different colors based on their weights. 
Light blue = light, yellow = average, Red = Heavy. 

Pressing "d" while on the simulation will display different mouse/keyboard interaction you can play around with.
Right clicking on the simulation will change the force between particle from Attraction to Repulsion (Interesting results
can be derived from this, my favorite is letting particles clump together and then "exploding them" by changing the force to
repulsion).

Pressing p enables/disables repulsion which is what enables or disables CCD (Continuous Collision Detection). My implementation
of CCD is not optimal as it is a brute force approach that doesn't take advantage of the QuadTree structures I'm using in
my Barnes-Hut implementation. That is why performance differences between Barnes-Hut and Naive Brute Force is better identified
by disabling CCD and letting particles fly around with no contact.

If you want to implement a Robust and better CCD, feel free to do so as that was my next step in this project.