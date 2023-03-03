/*
 * Body.cpp
 *
 *  Created on: Dec 15, 2014
 *  Updated on: Aug 23, 2022
 *      Author: ahueck
 */

#include "Body.h"
#include "Logger.h"
#include <math.h> 
namespace practical::nbody { 

Body::Body(const Vec2f& pos, const Vec2f& v, const float mass) : m_pos(pos), m_v(v), m_a(), m_mass(mass) {
}

void Body::update(float dt) {
  // DONE: update Body's speed "v" and position "pos" based on acceleration and timestep dt
  m_v = m_v + m_a * dt;
  m_pos = m_pos + m_v * dt;
} 

void Body::applyForces(const Body& other) {
  // DONE: update Body's acceleration based on "other"
  static const float eps = 0.0001f;
  static const float G   = 6.673e-11f;
  Vec2f delta = other.m_pos - m_pos; 
  float r = sqrt(pow(delta[0], 2) + pow(delta[1], 2) + eps);
  float fac = 1.0 / pow(r, 3);
  Vec2f m_a_change = delta * (G * other.m_mass * fac);
  m_a = m_a + m_a_change;
}

void Body::resetAcceleration() {
  m_a.reset();
}

const Vec2f& Body::position() const {
  return m_pos;
}

const Vec2f& Body::velocity() const {
  return m_v;
}

float Body::mass() const {
  return m_mass;
}

Vec2f& Body::velocity() {
  return m_v;
}

float& Body::mass() {
  return m_mass;
}

float Body::distanceTo(const Body& other) const {
  return m_pos.distance(other.m_pos);
}

std::ostream& operator<<(std::ostream& os, const Body& body) {
  const Vec2f& pos = body.position();
  const Vec2f& v   = body.velocity();
  os << "|p=(" << pos[0] << ", " << pos[1] << ") v=(" << v[0] << ", " << v[1] << ") m=" << body.mass() << "|";
  return os;
}

// hinzugefügte Setter-Methoden
void Body::setMass(const float new_mass) {
  m_mass = new_mass;
}
void Body::setSpeed(const Vec2f& new_speed) {
  m_v = new_speed;
}
void Body::setAccel(const Vec2f& new_accel) {
  m_a = new_accel;
}
}  // namespace practical::nbody
