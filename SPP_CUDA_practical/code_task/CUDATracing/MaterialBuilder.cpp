#include "MaterialBuilder.h"

MaterialBuilder* MaterialBuilder::_instance = 0;

MaterialBuilder::MaterialBuilder() {
}

MaterialBuilder* MaterialBuilder::Instance() {
  if (_instance == 0) {
    _instance = new MaterialBuilder();
  }
  return _instance;
}

Material MaterialBuilder::getMaterial(int mat) {
  Material m;

  if (mat == COPPER) {
    Color a_copper = Color(0.19f, 0.07f, 0.02f);
    Color d_copper = Color(0.7f, 0.27f, 0.08f);
    Color s_copper = Color(0.25f, 0.13f, 0.08f);
    m              = Material(a_copper, d_copper, s_copper, 12.8f);
  }

  else if (mat == OBSIDIAN) {
    Color a_obsidian = Color(0.05375f, 0.05f, 0.06625f);
    Color d_obsidian = Color(0.18275f, 0.17f, 0.22525f);
    Color s_obsidian = Color(0.332741f, 0.328634f, 0.346435f);
    m                = Material(a_obsidian, d_obsidian, s_obsidian, 38.4f);
  }

  else if (mat == SILVER) {
    Color a_silver = Color(0.19225f, 0.19225f, 0.19225f);
    Color d_silver = Color(0.50754f, 0.50754f, 0.50754f);
    Color s_silver = Color(0.508273f, 0.508273f, 0.508273f);
    m              = Material(a_silver, d_silver, s_silver, 51.2f);
  }

  else if (mat == GOLD) {
    Color a_gold = Color(0.24725f, 0.1995f, 0.0745f);
    Color d_gold = Color(0.75164f, 0.60648f, 0.22648f);
    Color s_gold = Color(0.628281f, 0.555802f, 0.366065f);
    m            = Material(a_gold, d_gold, s_gold, 51.2f);
  }

  else if (mat == EMERALD) {
    Color a_gold = Color(0.0215f, 0.1745f, 0.0215f);
    Color d_gold = Color(0.07568f, 0.61424f, 0.07568f);
    Color s_gold = Color(0.633f, 0.727811f, 0.633f);
    m            = Material(a_gold, d_gold, s_gold, 153.0f);
  }

  else {
    Color a_blue = Color(0.0f, 0.0f, 1.0f);
    Color d_blue = Color(0.0f, 0.0f, 1.0f);
    Color s_blue = Color(0.0f, 0.0f, 1.0f);
    m            = Material(a_blue, d_blue, s_blue, 100.0f);
  }

  return m;
}