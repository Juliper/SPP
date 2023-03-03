#ifndef LIST_H
#define LIST_H

typedef struct
{
	int item_count;
	int current_max_size;
	char growable;

	void **items;
	char **names;	
} lista;

void list_make(lista *listo, int size, char growable);
int list_add_item(lista *listo, void *item, char *name);
char* list_print_items(lista *listo);
void* list_get_name(lista *listo, char *name);
void* list_get_index(lista *listo, int indx);
void* list_get_item(lista *listo, void *item_to_find);
int list_find(lista *listo, char *name_to_find);
void list_delete_index(lista *listo, int indx);
void list_delete_name(lista *listo, char *name);
void list_delete_item(lista *listo, void *item);
void list_delete_all(lista *listo);
void list_print_list(lista *listo);
void list_free(lista *listo);

void test_list();
#endif