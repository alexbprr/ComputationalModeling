#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

float tfinal = 50;
float dt = 0.001;
float dy = 0.1;
float dx = 0.1; 

float xdim = 8;
float ydim = 6; 

int xsize;
int ysize;
int popsize = 6; 

float sn = 0.1;
float kn = 0.5;
float lnb = 0.5;
float rb = 0.05;
float mn = 0.2;
float Dn = 0.05;
float Db = 0.01;
float Xn = 0.05; 

float verifyDensity(float value){
    if (value < 0.) return 0.;
    else return value;
}

int getIndex(int x, int y, int pindex, int t){
    return pindex + x*popsize + y*(xsize*popsize) + t*(xsize*ysize*popsize);    
}

float chemotaxis(float *u, int x, int y, int pindex, int chindex, int t){
    float resx = 0, resy = 0, flux_left = 0, flux_right = 0;
    float p = u[getIndex(x,y,pindex,t)];
    float ch = u[getIndex(x,y,chindex,t)];
    if (x > 0){
        float pm1 = u[getIndex(x-1,y,pindex,t)];
        float chm1 = u[getIndex(x-1,y,chindex,t)];
        if(ch -  chm1 > 0)        
            flux_left = -((ch - chm1) * pm1)/dx;
        else        
            flux_left = -((ch - chm1) * p)/dx;        
    }
    if(x < xsize-1){
        float pp1 = u[getIndex(x+1,y,pindex,t)];
        float chp1 = u[getIndex(x+1,y,chindex,t)];
        if(chp1 - ch > 0)        
            flux_right = (chp1 - ch)*p /dx;         
        else        
            flux_right = (chp1 - ch)*pp1 /dx;        
    }
    resx = (flux_left + flux_right)/dx;

    flux_left = 0, flux_right = 0;
    if (y > 0) {
        float pm1 = u[getIndex(x,y-1,pindex,t)];
        float chm1 = u[getIndex(x,y-1,chindex,t)];
        if(ch - chm1 > 0)        
            flux_left = -(ch - chm1)*pm1 /dy;         
        else
            flux_left = -(ch - chm1)*p /dy;
    }
    if(y < ysize-1) {
        float pp1 = u[getIndex(x,y+1,pindex,t)];
        float chp1 = u[getIndex(x,y+1,chindex,t)];
        if(chp1 - ch > 0)        
            flux_right = (chp1 - ch)*p /dy; 
        else 
            flux_right = (chp1 - ch)*pp1 /dy; 
    }
    resy = (flux_left + flux_right)/dy;
    return resx + resy; 
}

float diffusionX(float *u, int x, int y, int pindex, int t){
    float ux = u[getIndex(x,y,pindex,t)], uxp1, uxm1;
    if (x == 0) {
        uxp1 = u[getIndex(x+1,y,pindex,t)];
        return (uxp1 - ux) /(dx*dx);
    }
    else if (x == (xsize - 1)) {
        uxm1 = u[getIndex(x-1,y,pindex,t)];
        return (uxm1 - ux)/(dx*dx);
    }
    else {   
        uxp1 = u[getIndex(x+1,y,pindex,t)];
        uxm1 = u[getIndex(x-1,y,pindex,t)];
        return (uxp1 - 2*ux + uxm1)/(dx*dx);
    }
}

float diffusionY(float *u, int x, int y, int pindex, int t){
    float uy = u[getIndex(x,y,pindex,t)], uyp1, uym1;
    if (y == 0) {
        uyp1 = u[getIndex(x,y+1,pindex,t)]; 
        return (uyp1 - uy) /(dy*dy);
    }
    else if (y == (ysize - 1)) {
        uym1 = u[getIndex(x,y-1,pindex,t)];
        return (uym1 - uy)/(dy*dy);
    }
    else {    
        uyp1 = u[getIndex(x,y+1,pindex,t)]; 
        uym1 = u[getIndex(x,y-1,pindex,t)];
        return (uyp1 - 2*uy + uym1)/(dy*dy);
    }
}

float diffusion(float* u, int x, int y, int pindex, int t){
    if(xsize == 1)
        return diffusionY(u,x,y,pindex,t);
    else if (ysize == 1)
        return diffusionX(u,x,y,pindex,t);
    else
        return diffusionX(u,x,y,pindex,t) + diffusionY(u,x,y,pindex,t);
}

void saveResults(FILE* arq, float* u, int p, int t){
    for (int y = ysize - 1; y >= 0; y--) {            
        for (int x = 0; x < xsize - 1; x++)
            fprintf(arq, "%lf , ", u[getIndex(x,y,p,t)]);
        fprintf(arq, "%lf ", u[getIndex(xsize-1,y,p,t)]);
        fprintf(arq,"\n");
    }   
    fprintf(arq,"\n"); 
}

typedef float * (*functiontype)(float, float, float, float, float, float);

float* ode(float v, float i, float e, float t, float b, float a){
    float *dudt = (float*) malloc(popsize*sizeof(float));
    //Parâmetros do modelo 
    dudt[0] = 0;
    dudt[1] = 0;
    dudt[2] = 0;
    dudt[3] = 0;
    dudt[4] = 0;
    dudt[5] = 0;
    return dudt;
}

float* rungeKutta(float *y0, float t, float h, functiontype f){
    float *y = y0;    
    float *k1 = f(y[0], y[1], y[2], y[3], y[4], y[5]);
    
    float *k2 = f(y[0] + (h*0.5)*k1[0], y[1] + (h*0.5)*k1[1], y[2] + (h*0.5)*k1[2], 
        y[3] + (h*0.5)*k1[3], y[4] + (h*0.5)*k1[4], y[5] + (h*0.5)*k1[5]);

    float *k3 = f(y[0] + (h*0.5)*k2[0], y[1] + (h*0.5)*k2[1], y[2] + (h*0.5)*k2[2], 
        y[3] + (h*0.5)*k2[3], y[4] + (h*0.5)*k2[4], y[5] + (h*0.5)*k2[5]);
        
    float *k4 = f(y[0] + h*k3[0], y[1] + h*k3[1], y[2] + h*k3[2], y[3] + h*k3[3], 
        y[4] + h*k3[4], y[5] + h*k3[5]); 
    
    for (int i = 0; i < popsize; i++)
        y[i] = y[i] + (h/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    return y;
}

void solveReactions(float *u, int x, int y, float t, int tnext, int tprev){        
    float *Y = (float*) malloc(popsize*sizeof(float));
    Y[0] = u[getIndex(x,y,0,tprev)];
    Y[1] = u[getIndex(x,y,1,tprev)];
    
    Y = rungeKutta(Y, t, dt, &(ode));
    
    u[getIndex(x,y,0,tnext)] = Y[0];
    u[getIndex(x,y,1,tnext)] = Y[1]; 
    
    free(Y);
}

void solveDiffusion(float *u, float *utemp, int x, int y, float t, int tnext){
    utemp[getIndex(x,y,0,tnext)] = u[getIndex(x,y,0,tnext)] + (Db*diffusion(u,x,y,0,tnext))*dt;
    utemp[getIndex(x,y,1,tnext)] = u[getIndex(x,y,1,tnext)] + (Dn*diffusion(u,x,y,1,tnext))*dt;
}

void solveChemotaxis(float *u, float *utemp, int x, int y, float t, int tnext){
    utemp[getIndex(x,y,1,tnext)] = u[getIndex(x,y,1,tnext)] - (Xn*chemotaxis(u,x,y,1,0,tnext))*dt;
}

void copyArray(float *u, float *utemp, int tnext){
    #pragma omp parallel for num_threads(6) 
    for (int y = 0; y < ysize; y++) {
        for (int x = 0; x < xsize; x++) {
            //Para cada populacao, copiar utemp para u
            for (int i = 0; i < popsize; i++)
                u[getIndex(x,y,i,tnext)] = utemp[getIndex(x,y,i,tnext)];
        }
    }
}

int main(){
    int tprev, tnext;         
    int nsteps = (int)((tfinal+dt)/dt);
    printf("T final: %.2f\n", tfinal);
    printf("steps: %d\n", nsteps);   
    int interval = nsteps/10;
    printf("interval: %d\n", interval); 

    xsize = (xdim/dx);
    printf("xsize: %d\n", xsize);    
    ysize = (ydim/dy);
    printf("ysize: %d\n", ysize);
    printf("domain size: %d\n", ysize*xsize);
        
    FILE *xFile = fopen("results/x.csv", "w"); 
    for (float x = 0; x < xdim; x += dx)
        fprintf(xFile, "%.3lf\n", x);
    fclose(xFile);

    FILE *yFile = fopen("results/y.csv", "w");
    for (float y = 0; y < ydim; y += dy)
        fprintf(yFile, "%.3lf\n", y);
    fclose(yFile);

    FILE *tFile = fopen("results/t.csv", "w");
    for (float t = 0; t < tfinal+dt; t += dt)
        fprintf(tFile, "%.6lf\n", t);
    fclose(tFile);

    int totalsize = (popsize)*ysize*xsize*2;
    printf("total size: %d\n", totalsize);
    float *u = (float*) malloc(totalsize*sizeof(float));
    float *utemp = (float*) malloc(totalsize*sizeof(float));

    for (int i = 0; i < totalsize; i++) {
        u[i] = 0;  
        utemp[i] = 0;
    } 

    u[getIndex(xsize/2,ysize/2,0,0)] = 100; 
    u[getIndex(xsize/2-1,ysize/2,0,0)] = 100; 
    u[getIndex(xsize/2+1,ysize/2,0,0)] = 100; 
    u[getIndex(xsize/2-5,ysize/2+5,0,0)] = 80;
    u[getIndex(10,10,0,0)] = 10; 
    u[getIndex(xsize-10,ysize-10,0,0)] = 50;
    u[getIndex(xsize-10,10,0,0)] = 30;

    tprev = tnext = 0;

    char bname[100], nname[100];
    sprintf(bname, "results/bacteria%d.csv", 0);
    sprintf(nname, "results/neutrophil%d.csv", 0);
    FILE *bFile = fopen(bname, "w");
    FILE *nFile = fopen(nname, "w");
    saveResults(bFile,u,0,0);
    saveResults(nFile,u,1,0);
    fclose(bFile);
    fclose(nFile);

    float t = 0;
    for(int step = 1; step <= nsteps; ++step) {

        if(step % 2 == 0){
            tnext = 0;
            tprev = 1;
        }
        else {
            tnext = 1;
            tprev = 0;
        }

        // Solving the PDEs with operator splitting: First solve reactions then the processes of movement                 
        #pragma omp parallel for num_threads(6) 
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {             
                solveReactions(u,x,y,t,tnext,tprev);                
            }
        }

        #pragma omp parallel for num_threads(6) 
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {             
                solveDiffusion(u,utemp,x,y,t,tnext);                
            }
        }

        copyArray(u,utemp,tnext);

        #pragma omp parallel for num_threads(6) 
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {             
                solveChemotaxis(u,utemp,x,y,t,tnext);                 
            }
        }

        copyArray(u,utemp,tnext);

        /*//Test 
        float totalB = 0;
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {             
                totalB += u[getIndex(x,y,0,tnext)];                
            }
        }
        printf("total de bacterias: %f\n", totalB);*/

        if (step % interval == 0) {
            char bname[100], nname[100];
            sprintf(bname, "results/bacteria%d.csv", step);
            sprintf(nname, "results/neutrophil%d.csv", step);
            FILE *bFile = fopen(bname, "w");
            FILE *nFile = fopen(nname, "w");
            saveResults(bFile,u,0,tnext);
            saveResults(nFile,u,1,tnext);
            fclose(bFile);
            fclose(nFile);
        }

        t += dt;  
    }

    free(u);
    free(utemp);
    return 0; 
}