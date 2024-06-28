
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER
#include <stdbool.h>
#include <pthread.h>
#else
#define bool int
#define true 1
#define false 0
#define strdup _strdup
#endif
#include <signal.h>

#include "ngspice_interface.h"

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


void start()
{
    int ret = 0;
    ret = ngSpice_Init(ng_getchar, ng_getstat, ng_exit,  ng_data, ng_initdata, ng_thread_runs, NULL);
    return;
}

void command(char* input)
{
    int ret = 0;
    ret = ngSpice_Command(input);
    return;
}

char** get_all_plots(){
    
    char** ret = ngSpice_AllPlots();
    int size = sizeof(ret);
    return ret;
}

pvector_info get_vector_info(char* vecname){
    pvector_info ret = ngGet_Vec_Info(vecname);
    return ret;
}

void circ(char** input){
    int ret = ngSpice_Circ(input);
    return;
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

    if(quitexit) {
        printf("DNote: Returned form quit with exit status %d\n", exitstatus);
        exit(exitstatus);
    }
    if(immediate) {
        printf("DNote: Unloading ngspice inmmediately is not possible\n");
        printf("DNote: Can we recover?\n");
    }

    else {
        printf("DNote: Unloading ngspice is not possible\n");
        printf("DNote: Can we recover? Send 'quit' command to ngspice.\n");
        errorflag = true;
        ngSpice_Command("quit 5");
//        raise(SIGINT);
    }

    return exitstatus;
}
