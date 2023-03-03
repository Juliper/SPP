#ifndef OBJPARSER_H
#define OBJPARSER_H

#include "list.h"

#define OBJ_FILENAME_LENGTH 500
#define MATERIAL_NAME_SIZE 255
#define OBJ_LINE_SIZE 500
#define MAX_VERTEX_COUNT 4 //can only handle quads or triangles

struct obj_face
{
	int vertex_index[MAX_VERTEX_COUNT];
	int normal_index[MAX_VERTEX_COUNT];
	int texture_index[MAX_VERTEX_COUNT];
	int vertex_count;
	int material_index;
};

struct obj_sphere
{
	int pos_index;
	int up_normal_index;
	int equator_normal_index;
	int texture_index[MAX_VERTEX_COUNT];
	int material_index;
};

struct obj_plane
{
	int pos_index;
	int normal_index;
	int rotation_normal_index;
	int texture_index[MAX_VERTEX_COUNT];
	int material_index;
};

struct obj_vector
{
	double e[3];
};

struct obj_material
{
	char name[MATERIAL_NAME_SIZE];
	char texture_filename[OBJ_FILENAME_LENGTH];
	double amb[3];
	double diff[3];
	double spec[3];
	double reflect;
	double refract;
	double trans;
	double shiny;
	double glossy;
	double refract_index;
};

struct obj_camera
{
	int camera_pos_index;
	int camera_look_point_index;
	int camera_up_norm_index;
};

struct obj_light_point
{
	int pos_index;
	int material_index;
};

struct obj_light_disc
{
	int pos_index;
	int normal_index;
	int material_index;
};

struct obj_light_quad
{
	int vertex_index[MAX_VERTEX_COUNT];
	int material_index;
};

struct obj_growable_scene_data
{
//	vector extreme_dimensions[2];
	char scene_filename[OBJ_FILENAME_LENGTH];
	char material_filename[OBJ_FILENAME_LENGTH];
	
	lista vertex_list;
	lista vertex_normal_list;
	lista vertex_texture_list;
	
	lista face_list;
	lista sphere_list;
	lista plane_list;
	
	lista light_point_list;
	lista light_quad_list;
	lista light_disc_list;
	
	lista material_list;
	
	obj_camera *camera;

	float maxCoordVal;
};

struct obj_scene_data
{
	obj_vector **vertex_list;
	obj_vector **vertex_normal_list;
	obj_vector **vertex_texture_list;
	obj_vector **face_normal_list;	
	obj_face **face_list;
	obj_sphere **sphere_list;
	obj_plane **plane_list;
	
	obj_light_point **light_point_list;
	obj_light_quad **light_quad_list;
	obj_light_disc **light_disc_list;
	
	obj_material **material_list;
	
	int vertex_count;
	int vertex_normal_count;
	int vertex_texture_count;

	int face_count;
	int sphere_count;
	int plane_count;

	int light_point_count;
	int light_quad_count;
	int light_disc_count;

	int material_count;

	obj_camera *camera;
};







int parse_obj_scene(obj_scene_data *data_out, char *filename);
void delete_obj_data(obj_scene_data *data_out);
void obj_free_half_list(lista *listo);
int obj_convert_to_list_index(int current_max, int index);
void obj_convert_to_list_index_v(int current_max, int *indices);
void obj_set_material_defaults(obj_material *mtl);
int obj_parse_vertex_index(int *vertex_index, int *texture_index, int *normal_index);
obj_face* obj_parse_face(obj_growable_scene_data *scene);
obj_sphere* obj_parse_sphere(obj_growable_scene_data *scene);
obj_plane* obj_parse_plane(obj_growable_scene_data *scene);
obj_light_point* obj_parse_light_point(obj_growable_scene_data *scene);
obj_light_quad* obj_parse_light_quad(obj_growable_scene_data *scene);
obj_light_disc* obj_parse_light_disc(obj_growable_scene_data *scene);
obj_vector* obj_parse_vector();
void obj_parse_camera(obj_growable_scene_data *scene, obj_camera *camera);
int obj_parse_mtl_file(char *filename, lista *material_list);
int obj_parse_obj_file(obj_growable_scene_data *growable_data, char *filename);
void obj_init_temp_storage(obj_growable_scene_data *growable_data);
void obj_free_temp_storage(obj_growable_scene_data *growable_data);
void delete_obj_data(obj_scene_data *data_out);
void obj_copy_to_out_storage(obj_scene_data *data_out, obj_growable_scene_data *growable_data);
#endif 