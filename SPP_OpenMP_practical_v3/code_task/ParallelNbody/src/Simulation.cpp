/*
 * Simulation.cpp
 *
 *  Created on: Dec 15, 2014
 *  Updated on: Aug 23, 2022
 *      Author: ahueck
 */

#include "Simulation.h"
#include "Logger.h"
#include "Vec2.h"
#include <omp.h>

namespace practical::nbody {

Simulation::Simulation(SimulationData data) : m_sim_data(std::move(data)) {
}

void Simulation::run() {
  for (unsigned int step = 0; step < m_sim_data.num_steps; ++step) {
    nextTimestep();
  }
}

void Simulation::nextTimestep() {
   // DONE: implement N-body calculation according to the exercise sheet.

  // Kollisionen erkennen und verarbeiten
  handleCollisions();

  // Beschleunigungen der Körper aktualisieren (basierend
  // auf ihren gegenseitigen Anziehungskräften)
  #pragma omp parallel for
  for (unsigned int i = 0; i < m_sim_data.bodies.size(); i++) {
    Body* body_i = &m_sim_data.bodies[i];
    body_i->setAccel({0.0, 0.0});
    for (unsigned int j = 0; j < m_sim_data.bodies.size(); j++) {
      
      Body* body_j = &m_sim_data.bodies[j];
      if (i != j && (int) body_i->mass() != 0 && (int) body_j->mass() != 0) {
        body_i->applyForces(body_j[0]); 
      }
    }
  }

  // Geschwindigkeit und Position der Körper  
  // auf Grundlage der neuen Beschleunigung aktualisieren.
  for (unsigned int i = 0; i < m_sim_data.bodies.size(); i++) {
    Body* body = &m_sim_data.bodies[i];
    if (body->mass() != 0.0) {
      body->update(m_sim_data.dt);
    }
   } 
}

void Simulation::handleCollisions() {
  // DONE: handle collision per timestep (before updating forces)  
  std::vector<std::vector<Body*>> collision_list;

  // Alle Paare kollidierender Körper im Vektor collision_list speichern
  #pragma omp parallel for 
  for (unsigned int i = 0; i < m_sim_data.bodies.size(); i++) {
    for (unsigned int j = 0; j < m_sim_data.bodies.size(); j++) {
      Body* body_i = &m_sim_data.bodies[i];
      Body* body_j = &m_sim_data.bodies[j];
      if (i < j && body_i->mass() != 0.0 && body_j->mass() != 0.0) {
        // Kollisionskontrolle
        if (body_i->distanceTo(body_j[0]) < m_sim_data.distance) {
           const std::vector<Body*> &colliding_pair = {body_i, body_j};
           #pragma omp critical 
           {
            collision_list.push_back(colliding_pair);
           }
        }
      }
    }
  }

  // Alle Kollisionen verarbeiten
  for (std::vector<Body*>& colliding_pair: collision_list) {
    Body* body_i = colliding_pair[0];
    Body* body_j = colliding_pair[1];
    if(body_i->mass() != 0 && body_j->mass() != 0) {
      float mass_ges = body_i->mass() + body_j->mass();
      body_i->setSpeed((body_i->velocity() * body_i->mass() + body_j->velocity() * body_j->mass()) / mass_ges);
      body_i->setMass(mass_ges);
      body_j->setMass(0.0);
    }
  }
  
}

const SimulationData& Simulation::simulationState() const {
  return m_sim_data;
}

Simulation::~Simulation() = default;

}  // namespace practical::nbody