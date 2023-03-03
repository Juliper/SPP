#include "rtmLoader.h"

#include "string_util.h"

void path_from_filename(char* path, int pathlen, const char* filename) {
  const char* ptr     = strrchr(filename, '/');
  const char* ptr_bsl = strrchr(filename, '\\');
  if (ptr_bsl > ptr)
    ptr = ptr_bsl;

  pathlen--;
  int l = (int)(ptr - filename);
  if (l > pathlen)
    l = pathlen;
  memcpy(path, filename, l);
  path[l] = 0;
}

rtmLoader::rtmLoader() {
}

Mesh rtmLoader::loadRtm(const std::string& filename, Transformation* trans) {
  printf("Load %s\n", filename.c_str());

  const int maxLineLength = 512;
  char line[maxLineLength];
  int lineNumber = 0, nTabs, l;
  int err        = 0;
  FILE* F;
  char entryname[256];
  int i1, vtxIdx[3];

  const int path_len = 1024;
  char path[path_len];
  path_from_filename(path, path_len, filename.c_str());

  Mesh mesh = Mesh();

  if ((F = fopen(filename.c_str(), "rb")) == NULL) {
    // throw "scene::load_rtm(): can not open '" + filename + "'";
  }

  while (!feof(F)) {
    if ((l = read_line(F, line, maxLineLength, &nTabs)) < 0) {
      break;
    }
    lineNumber++;
    if (empty_line(line))
      continue;

    if (sscanf(line, "%s ", entryname) != 1)
      continue;

    if (!strcmpnc(entryname, "number_of_faces")) {  // number_of_faces
      if ((sscanf(line, "%s %d", entryname, &i1) != 2) || (i1 < 0)) {
        err = -1;
        break;
      }
      mesh.setFaceCount(i1);
    }

    else if (!strcmpnc(entryname, "number_of_vertices")) {  // number_of_vertices
      if ((sscanf(line, "%s %d", entryname, &i1) != 2) || (i1 < 0)) {
        err = -1;
        break;
      }
      mesh.setVtxCount(i1);
    }

    else if (!strcmpnc(entryname, "face")) {  // face
      if (sscanf(line, "%s %d vtx %d %d %d", entryname, &i1, vtxIdx + 0, vtxIdx + 1, vtxIdx + 2) != 5) {
        err = -1;
        break;
      }
      mesh.addFace(i1, Face(vtxIdx[0], vtxIdx[1], vtxIdx[2]));
    }

    else if (!strcmpnc(entryname, "vertex")) {  // vertex
      float x, y, z, nx, ny, nz, u, v;
      if ((sscanf(line, "%s %d pos %f %f %f normal %f %f %f uv %f %f", entryname, &i1, &x, &y, &z, &nx, &ny, &nz, &u,
                  &v) != 10)) {
        err = -1;
        break;
      }

      Point p     = Point(x, y, z);
      Direction n = Direction(nx, ny, nz);
      p.n         = n;
      p.u         = u;
      p.v         = v;

      mesh.addPoint(i1, trans->transformPoint(p));
    }
  }

  if (err) {
    // error detected:
    char buf[1024];
    // snprintf(buf, "scene::load_rtm(): Invalid entry '%s' at line %d\n", entryname, lineNumber);
    throw std::string(buf);
  }

  fclose(F);
  return mesh;
}
