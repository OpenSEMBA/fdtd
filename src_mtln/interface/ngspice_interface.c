
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ngspice/bool.h"

#ifndef _MSC_VER
#include <pthread.h>
#else
#undef bool
#define bool int
#define true 1
#define false 0
#define strdup _strdup
#endif
#include <signal.h>
#include "ngspice_interface.h"
#define ngspice_COMPLEX_H
#include "ngspice/fteext.h"
#include "ngspice/ftedefs.h"
#ifndef _WIN32
#include <unistd.h>
#endif

int INPtypelook(char *);
int load_opus(const char *);

#ifndef _MSC_VER
#ifndef __cplusplus
typedef int bool;
#endif
typedef int BOOL;
#define true 1
#define false 0
#endif

bool no_bg = true;
static bool errorflag = false;

int
ng_getchar(char* outputreturn, int ident, void* userdata);

int
ng_getstat(char* outputreturn, int ident, void* userdata);

int
ng_exit(int, bool, bool, int ident, void*);

int
ng_thread_runs(bool noruns, int ident, void* userdata);

int
ng_initdata(pvecinfoall intdata, int ident, void* userdata);

int
ng_data(pvecvaluesall vdata, int numvecs, int ident, void* userdata);


/* Folder of the running executable, with trailing separator; empty on failure. */
static void executable_folder(char* folder, size_t size)
{
    char* sep;
    folder[0] = '\0';
#ifdef _WIN32
    char* path = NULL;
    if (_get_pgmptr(&path) != 0 || path == NULL) return;
    strncpy(folder, path, size - 1);
    folder[size - 1] = '\0';
    sep = strrchr(folder, '\\');
    if (strrchr(folder, '/') > sep) sep = strrchr(folder, '/');
#else
    ssize_t len = readlink("/proc/self/exe", folder, size - 1);
    if (len <= 0) { folder[0] = '\0'; return; }
    folder[len] = '\0';
    sep = strrchr(folder, '/');
#endif
    if (sep) sep[1] = '\0'; else folder[0] = '\0';
}

/* XSPICE code models (e.g. 'filesource' used by cable-end generators) are
   shipped next to the executable. They are loaded only if no spinit has
   already loaded them, so a user spinit keeps precedence. */
static void load_codemodels_next_to_executable(void)
{
    static const char* models[] = {"spice2poly.cm", "analog.cm", "digital.cm",
                                   "xtradev.cm", "xtraevt.cm", "table.cm"};
    char folder[4096], path[4200];
    size_t i;
    FILE* f;

    if (INPtypelook("filesource") >= 0) return;
    executable_folder(folder, sizeof(folder));
    if (folder[0] == '\0') return;
    for (i = 0; i < sizeof(models) / sizeof(models[0]); i++) {
        snprintf(path, sizeof(path), "%s%s", folder, models[i]);
        f = fopen(path, "rb");
        if (!f) continue;
        fclose(f);
        load_opus(path);
    }
}

void start()
{
    int ret = 0;
    errorflag = false;
    ret = ngSpice_Init(ng_getchar, ng_getstat, ng_exit,  ng_data, ng_initdata, ng_thread_runs, NULL);
    load_codemodels_next_to_executable();
    return;
}

void command(char* input)
{
    int ret = 0;
    if (errorflag) {
        return;
    }
    ret = ngSpice_Command(input);
    return;
}

char** get_all_plots(){
    
    char** ret = ngSpice_AllPlots();
    int size = sizeof(ret);
    return ret;
}

pvector_info get_vector_info(char* vecname){
    if (errorflag) {
        return NULL;
    }
    pvector_info ret = ngGet_Vec_Info(vecname);
    return ret;
}

/* Returns 0 if the circuit was parsed, non-zero otherwise.
   ngSpice_Circ returns 0 even when parsing fails, so the parsed circuit is checked. */
int circ(char** input){
    if (errorflag) {
        return 1;
    }
    if (ngSpice_Circ(input) != 0) return 1;
    return (ft_curckt && ft_curckt->ci_ckt) ? 0 : 1;
}

int has_error(void)
{
    return errorflag ? 1 : 0;
}

/* Callback function called from bg thread in ngspice to transfer
any string created by printf or puts. Output to stdout in ngspice is
preceded by token stdout, same with stderr.*/
int
ng_getchar(char* outputreturn, int ident, void* userdata)
{
    // printf("%s\n", outputreturn);
    /* setting a flag if an error message occurred */
    return 0;
}

int
ng_getstat(char* outputreturn, int ident, void* userdata)
{
    // printf("%s\n", outputreturn);
    return 0;
}

int
ng_thread_runs(bool noruns, int ident, void* userdata)
{
    no_bg = noruns;
    return 0;
}

/* Callback function called from bg thread in ngspice once per accepted data point */
int
ng_data(pvecvaluesall vdata, int numvecs, int ident, void* userdata)
{
    return 0;
}


int
ng_initdata(pvecinfoall intdata, int ident, void* userdata)
{
    return 0;
}


/* Callback function called from bg thread in ngspice if fcn controlled_exit()
   is hit. Do not exit, but unload ngspice. */
int
ng_exit(int exitstatus, bool immediate, bool quitexit, int ident, void* userdata)
{
    errorflag = true;

    if(quitexit) {
        printf("DNote: ngspice requested quit with exit status %d\n", exitstatus);
        return exitstatus;
    }
    if(immediate) {
        printf("DNote: Unloading ngspice inmmediately is not possible\n");
        printf("DNote: Can we recover?\n");
    }

    else {
        printf("DNote: Unloading ngspice is not possible\n");
        printf("DNote: Can we recover? Send 'quit' command to ngspice.\n");
    }

    return exitstatus;
}


