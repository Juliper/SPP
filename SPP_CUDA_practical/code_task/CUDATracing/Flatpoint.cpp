#include "Flatpoint.h"

#include "Constants.h"
#include "Log.h"

FlatPoint::FlatPoint(Plane* pl, Ray* r) {
  // TEST
  /*r->ray7 = 0.11f;
  r->ray8 = -0.22f;
  r->ray9 = 0.33f;
  r->ray11 = -0.44f;
  r->ray12 = 0.55;
  r->ray14 = -0.66f;*/

  this->fp9  = -pl->pl5 * r->ray11 + pl->pl4 * r->ray12 - pl->pl3 * r->ray14;  // e1^einf
  this->fp12 = pl->pl5 * r->ray8 - pl->pl4 * r->ray9 + pl->pl2 * r->ray14;     // e2^einf
  this->fp14 = -pl->pl5 * r->ray7 + pl->pl3 * r->ray9 - pl->pl2 * r->ray12;    // e3^einf
  this->fp16 = pl->pl4 * r->ray7 + pl->pl2 * r->ray11 - pl->pl3 * r->ray8;     // einf^e0

  this->fp9  = 1.0f / (-fp16) * fp9;   // e1
  this->fp12 = 1.0f / (-fp16) * fp12;  // e2
  this->fp14 = 1.0f / (-fp16) * fp14;  // e3
  this->fp16 = 1.0f;                   // e0
}
