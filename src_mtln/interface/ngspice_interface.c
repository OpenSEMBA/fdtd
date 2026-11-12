
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
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

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
/* Only the sharedspice API is exported by ngspice.dll on Windows, so ngspice
   state is probed by matching the text a command prints. */
static const char* probe_text = NULL;
static bool probe_hit = false;

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

/* Runs an ngspice command and reports whether its output contained text. */
static bool command_prints(const char* cmd, const char* text)
{
    char buffer[4300];
    snprintf(buffer, sizeof(buffer), "%s", cmd);
    probe_text = text;
    probe_hit = false;
    ngSpice_Command(buffer);
    probe_text = NULL;
    return probe_hit;
}

/* XSPICE code models (e.g. 'filesource' used by cable-end generators) are
   shipped next to the executable. They are loaded only if no spinit has
   already loaded them, so a user spinit keeps precedence. */
static void load_codemodels_next_to_executable(void)
{
    static const char* models[] = {"spice2poly.cm", "analog.cm", "digital.cm",
                                   "xtradev.cm", "xtraevt.cm", "table.cm"};
    char folder[4096], path[4200], cmd[4300];
    size_t i;
    FILE* f;

    if (!command_prints("devhelp filesource", "not found")) return;
    executable_folder(folder, sizeof(folder));
    if (folder[0] == '\0') return;
    for (i = 0; i < sizeof(models) / sizeof(models[0]); i++) {
        snprintf(path, sizeof(path), "%s%s", folder, models[i]);
        f = fopen(path, "rb");
        if (!f) continue;
        fclose(f);
#ifdef _WIN32
        /* ngspice splits command words on spaces; the 8.3 path has none. */
        {
            char short_path[4200];
            DWORD n = GetShortPathNameA(path, short_path, sizeof(short_path));
            if (n > 0 && n < sizeof(short_path)) strcpy(path, short_path);
        }
#endif
        snprintf(cmd, sizeof(cmd), "codemodel %s", path);
        ngSpice_Command(cmd);
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
    return command_prints("inventory", "no current circuit") ? 1 : 0;
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
    if (probe_text && strstr(outputreturn, probe_text)) probe_hit = true;
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


