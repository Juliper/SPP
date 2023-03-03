#include "MeshFileReader.h"

#include <exception>
#include <math.h>

MeshFileReader::MeshFileReader() {
}

Mesh MeshFileReader::readMeshFromFile(char* filename, Transformation* trans, int filetype) {
  // obj
  if (filetype == 1)
    return readObjFile(filename, trans);
  // rtm
  else if (filetype == 2)
    return readRtmFile(filename, trans);
  else
    return readRtmFile("files\\teapot.rtm", trans);
}

Mesh MeshFileReader::readObjFile(char* filename, Transformation* trans) {
  float MaxCoordVal  = 0.0f;
  objLoader* objLoad = new objLoader();
  objLoad->load(filename);

  Mesh mesh = Mesh();

  // Punkte kopieren
  int vtxcnt = objLoad->vertexCount;
  mesh.setVtxCount(vtxcnt);
  for (int i = 0; i < vtxcnt; i++) {
    Point p = Point((float)objLoad->vertexList[i]->e[0], (float)objLoad->vertexList[i]->e[1],
                    (float)objLoad->vertexList[i]->e[2]);

    if (abs(p.point2) > MaxCoordVal)
      MaxCoordVal = abs(p.point2);

    if (abs(p.point3) > MaxCoordVal)
      MaxCoordVal = abs(p.point3);

    if (abs(p.point2) > MaxCoordVal)
      MaxCoordVal = abs(p.point3);

    Direction n;
    if (objLoad->normalCount == objLoad->vertexCount) {
      n = Direction((float)objLoad->normalList[i]->e[0], (float)objLoad->normalList[i]->e[1],
                    (float)objLoad->normalList[i]->e[2]);
      n.normalize();
    }

    p.n = n;
    mesh.addPoint(i, p);
  }
  // Punkte normalisieren
  for (int i = 0; i < vtxcnt; i++) {
    mesh.points[i].point2 = mesh.points[i].point2 / MaxCoordVal;
    mesh.points[i].point3 = mesh.points[i].point3 / MaxCoordVal;
    mesh.points[i].point4 = mesh.points[i].point4 / MaxCoordVal;

    mesh.points[i] = trans->transformPoint(mesh.points[i]);
  }

  // Faces kopieren
  int facecnt = objLoad->faceCount;
  mesh.setFaceCount(facecnt);
  for (int i = 0; i < facecnt; i++) {
    mesh.addFace(i, Face(objLoad->faceList[i]->vertex_index[0], objLoad->faceList[i]->vertex_index[1],
                         objLoad->faceList[i]->vertex_index[2]));
  }

  if (objLoad->normalCount != objLoad->vertexCount)
    mesh.interpolateNormals();

  delete objLoad;
  return mesh;
}

Mesh MeshFileReader::readRtmFile(char* filename, Transformation* trans) {
  rtmLoader* reader = new rtmLoader();
  Mesh m            = reader->loadRtm(filename, trans);
  delete reader;
  return m;
}