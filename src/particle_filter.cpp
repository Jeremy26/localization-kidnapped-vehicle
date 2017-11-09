/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	num_particles = 100;
	cout<<"INIT"<<endl;
	
	default_random_engine gen;
	// Add random Gaussian noise to each particle.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	for (int i = 0; i<num_particles;i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;
		particles.push_back(particle);
		weights.push_back(1);

		cout<<"["<<particle.x<<" ;"<<particle.y<<" ;"<<particle.theta<<" ;"<<particle.weight<<" ]"<<endl;
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	cout<<"PREDICT"<<endl;
	
	for (int i=0;i<particles.size();i++)
	{
		double theta_temp = particles[i].theta;

		if (yaw_rate == 0){
			particles[i].x += velocity*delta_t*cos(theta_temp);
			particles[i].y += velocity*delta_t*sin(theta_temp);
		}
		else {
			particles[i].x += velocity*(sin(theta_temp + delta_t*yaw_rate) - sin(theta_temp))/yaw_rate;
			particles[i].y += velocity*(cos(theta_temp) - cos(theta_temp+yaw_rate*delta_t))/yaw_rate;
			particles[i].theta += yaw_rate*delta_t;
		}

		normal_distribution<double> noise_x(particles[i].x,std_pos[0]);
		normal_distribution<double> noise_y(particles[i].y,std_pos[1]);
		normal_distribution<double> noise_theta(particles[i].theta, std_pos[2]);
		
		particles[i].x += noise_x(gen);
		particles[i].y += noise_y(gen);
		particles[i].theta += noise_theta(gen);
		
		cout<<"["<<particles[i].x<<" ;"<<particles[i].y<<" ;"<<particles[i].theta<<" ;"<<particles[i].weight<<" ]"<<endl;
		
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//## FOR EACH LANDMARK, GO THROUGH AND FIND THE CLOSEST MEASUREMENT
	
	double min_dist = 1000000; // We set up a high initial distance

	// Take the closest distance between measurements and observations and define observation id to be map id
	for (int i=0; i<predicted.size();i++){
		for (int j = 0; j<observations.size();j++)
		{
			if (dist(observations[j].x, observations[j].y, predicted[i].x, predicted[i].y) < min_dist)
			observations[j].id = predicted[i].id;
		}
	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	cout<<"UPDATE"<<endl;
	
	for (int p=0; p<num_particles;p++){
		//############# STEP 1 : MATCH MAP LANDMARKS AND OBSERVATIONS COORDINATES  #############
		/*
		Our goal here is to arrive at step 2 : NEAREST NEIGHBOR
		That step will call the above DataAssociation function that takes in 2 LandmarkObs vectors.
		We take the map landmark observations and "cast" them into the LandmarkObs vector predictions we just created.

		Since the Map coordinates here are the reference coordinate, no need for any transformation
		*/
		vector<LandmarkObs> actual_landmarks; // Vector of LandmarksObs

		for(int i=0; i<map_landmarks.landmark_list.size();i++)
		{
			int landmark_id = map_landmarks.landmark_list[i].id_i;
			float landmark_x = map_landmarks.landmark_list[i].x_f;
			float landmark_y = map_landmarks.landmark_list[i].y_f;
			
			LandmarkObs land;
			land.x = landmark_x;
			land.y = landmark_y;
			land.id = landmark_id;
			actual_landmarks.push_back(land);
		}

		/*
		We now have one element ready for the dataAssociation function that takes in map landmarks and observations.
		Observations sent to that function will be subject to transformations. We do rotations and translations according to the course to do so.
		We iterate for each observation and do the transformation to have a vector of transformed observations.
		*/
		LandmarkObs obs, trans_obs; // Current observation and transformed observation
		vector<LandmarkObs> transformed_observations; // Vector for all transformed observations

		for (int i =0; i<observations.size();i++)
		{
			obs = observations[i];
			double theta_temp = particles[p].theta;
			trans_obs.x = particles[p].x + cos (theta_temp)*obs.x - sin(theta_temp)*obs.y;
			trans_obs.y = particles[p].y + sin(theta_temp)*obs.x - cos(theta_temp)*obs.y;
			transformed_observations.push_back(trans_obs);
		}

		// WE NOW HAVE 2 VECTOR OF LANDMARK OBS CONTAINING OBSERVATIONS AND ACTUAL MAP LANDMARKS
		//############# STEP 2 : NEAREST NEIGHBOR TECHNIQUE  #############
		
		dataAssociation(actual_landmarks,transformed_observations);

		// WE NOW HAVE OBSERVATIONS OF WHAT OUR LANDMARKS CAN BE		
		//############# STEP 3 : WEIGHT UPDATES USING MULTIVARIATE GAUSSIAN DISTRIBUTION  #############
		particles[p].weight = 1.0;
	  	
	  	vector<int> associations;
	  	vector<double> sense_x;
	  	vector<double> sense_y;
		
		for (int i =0; i<transformed_observations.size();i++)
		{
			double measured_x = transformed_observations[i].x;
			double measured_y = transformed_observations[i].y;
			double mu_x = map_landmarks.landmark_list[i].x_f;
			double mu_y = map_landmarks.landmark_list[i].y_f;
			double normalizer = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
			double exponent = pow(measured_x - mu_x,2)/(2*pow(std_landmark[0],2)) + pow(measured_y - mu_y,2)/(2*pow(std_landmark[1],2));
			double multiplier = normalizer*exp(-exponent);
			if (multiplier >0)
			{
				particles[p].weight *= multiplier;
			}
			associations.push_back(transformed_observations[i].id);
			sense_x.push_back(transformed_observations[i].x);
			sense_y.push_back(transformed_observations[i].y);
		}
		SetAssociations(particles[p], associations,sense_x,sense_y);
		
		weights[p] = particles[p].weight;
		cout<<"["<<particles[p].x<<" ;"<<particles[p].y<<" ;"<<particles[p].theta<<" ;"<<particles[p].weight<<" ]"<<endl;
		
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(),weights.end());

	vector<Particle> resample_particles;

	for(int i=0; i<num_particles;i++)
	{
		resample_particles.push_back(particles[distribution(gen)]);
	}
	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
