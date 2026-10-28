/*
 * planets_mpi.c  --  simple parallel N-body simulation of the solar system
 *
 *   reads  : initial.csv   (name,parent,gm_km3_s2,x,y,z,vx,vy,vz)
 *   writes : traj.csv      (t_days,name,x,y,z,vx,vy,vz)
 *
 *   run for 5 years, velocity Verlet, units AU / day / AU^3/day^2
 *
 * build : mpicc -O2 -o planets_mpi planets_mpi.c -lm
 * run   : mpirun -np 4 ./planets_mpi initial.csv traj.csv
 *
 * NOTE: with only ~40 bodies MPI is *slower* than a serial build.
 *       This is a learning exercise, not a speed-up.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define MAX_BODIES   128
#define MAX_NAME     64
#define AU_KM        149597870.7
#define DAY_S        86400.0
#define YEARS        5.0
#define DT           0.0001              /* days per step  */
#define SAVE_EVERY   10000             /* every 100 steps = 10 days */

typedef struct {
    char   name[MAX_NAME];
    double gm;                        /* AU^3 / day^2 */
    double x, y, z;                   /* AU           */
    double vx, vy, vz;                /* AU / day     */
} Body;

static Body   all[MAX_BODIES];        /* same list on every rank */
static double gm_all[MAX_BODIES];
static int    n = 0;

/* ------------------------------------------------------------------ */
/*  read initial.csv                                                  */
/* ------------------------------------------------------------------ */
static int read_initial(const char *fn)
{
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); return -1; }

    char line[1024];
    int  first = 1;

    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;

        if (first) { first = 0; continue; }               /* skip header */

        int commas = 0;
        for (char *p = line; *p; ++p) if (*p == ',') commas++;

        char   name[MAX_NAME], parent[MAX_NAME];
        double gm_km, x, y, z, vx, vy, vz;

        if (commas >= 9) {                                /* 10-column file */
            double gm_au_unused;
            if (sscanf(line, "%63[^,],%63[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       name, parent, &gm_km, &gm_au_unused,
                       &x, &y, &z, &vx, &vy, &vz) < 10) continue;
        } else {                                          /*  9-column file */
            if (sscanf(line, "%63[^,],%63[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       name, parent, &gm_km,
                       &x, &y, &z, &vx, &vy, &vz) < 9) continue;
        }

        Body *b = &all[n];
        strncpy(b->name, name, MAX_NAME - 1);
        b->gm = gm_km * (DAY_S * DAY_S) / (AU_KM * AU_KM * AU_KM);
        b->x = x;  b->y = y;  b->z = z;
        b->vx = vx; b->vy = vy; b->vz = vz;
        n++;
        if (n >= MAX_BODIES) break;
    }
    fclose(f);
    return n;
}

/* ------------------------------------------------------------------ */
/*  acceleration for the bodies owned by this rank                     */
/* ------------------------------------------------------------------ */
static void accel(const double *X, const double *Y, const double *Z,
                  const double *GM, int ntot,
                  int start, int count,
                  double *AX, double *AY, double *AZ)
{
    for (int i = 0; i < count; ++i) {
        int    I  = start + i;
        double ax = 0, ay = 0, az = 0;
        for (int J = 0; J < ntot; ++J) {
            if (J == I) continue;
            double dx = X[J] - X[I];
            double dy = Y[J] - Y[I];
            double dz = Z[J] - Z[I];
            
            // إضافة 1e-12 لمنع القسمة على صفر نهائياً
            double r2 = dx*dx + dy*dy + dz*dz + 1e-12; 
            
            double k  = GM[J] / (r2 * sqrt(r2));
            ax += k * dx;  ay += k * dy;  az += k * dz;
        }
        AX[i] = ax;  AY[i] = ay;  AZ[i] = az;
    }
}

/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const char *infile  = (argc > 1) ? argv[1] : "initial.csv";
    const char *outfile = (argc > 2) ? argv[2] : "traj.csv";

    if (rank == 0) read_initial(infile);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (n <= 0) { if (rank==0) fprintf(stderr,"no bodies\n"); MPI_Finalize(); return 1; }

    /* split bodies across ranks */
    int base = n / size, extra = n % size;
    int my_start = rank * base + (rank < extra ? rank : extra);
    int my_count = base + (rank < extra ? 1 : 0);

    for (int i = 0; i < n; ++i) gm_all[i] = all[i].gm;

    /* global state (every rank keeps a full copy) */
    double *X  = calloc(n, sizeof *X), *Y  = calloc(n, sizeof *Y), *Z  = calloc(n, sizeof *Z);
    double *VX = calloc(n, sizeof *VX), *VY = calloc(n, sizeof *VY), *VZ = calloc(n, sizeof *VZ);

    /* local work arrays */
    double *axl = calloc(my_count, sizeof *axl);
    double *ayl = calloc(my_count, sizeof *ayl);
    double *azl = calloc(my_count, sizeof *azl);
    double *vxh = calloc(my_count, sizeof *vxh);
    double *vyh = calloc(my_count, sizeof *vyh);
    double *vzh = calloc(my_count, sizeof *vzh);

    /* Allgatherv offsets */
    int *cnts = malloc(size * sizeof *cnts);
    int *offs = malloc(size * sizeof *offs);
    for (int r = 0; r < size; ++r) {
        cnts[r] = base + (r < extra ? 1 : 0);
        offs[r] = r * base + (r < extra ? r : extra);
    }

    /* initial state */
    for (int i = 0; i < n; ++i) {
        X[i]=all[i].x;  Y[i]=all[i].y;  Z[i]=all[i].z;
        VX[i]=all[i].vx; VY[i]=all[i].vy; VZ[i]=all[i].vz;
    }

    accel(X, Y, Z, gm_all, n, my_start, my_count, axl, ayl, azl);

    FILE *out = NULL;
    if (rank == 0) {
        out = fopen(outfile, "w");
        if (!out) { fprintf(stderr, "cannot open %s\n", outfile); MPI_Abort(MPI_COMM_WORLD,1); }
        fprintf(out, "t_days,name,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud\n");
    }

    int    n_steps = (int)(YEARS * 365.25 / DT + 0.5);
    double t       = 0.0;

    for (int step = 0; step <= n_steps; ++step) {

        /* ---------- write snapshot ---------------------------------- */
        if (step % SAVE_EVERY == 0) {
            /* sync velocities so rank 0 sees everybody */
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VX, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VY, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VZ, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);

            if (rank == 0) {
                for (int i = 0; i < n; ++i)
                    fprintf(out, "%.3f,%s,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                            t, all[i].name, X[i], Y[i], Z[i], VX[i], VY[i], VZ[i]);
            }
        }

        if (step == n_steps) break;

        /* ---------- velocity Verlet -------------------------------- */
        /* 1) half kick  v <- v + a * dt/2 */
        for (int i = 0; i < my_count; ++i) {
            vxh[i] = VX[my_start+i] + 0.5 * DT * axl[i];
            vyh[i] = VY[my_start+i] + 0.5 * DT * ayl[i];
            vzh[i] = VZ[my_start+i] + 0.5 * DT * azl[i];
        }

        /* 2) drift  x <- x + v_half * dt */
        for (int i = 0; i < my_count; ++i) {
            int I = my_start + i;
            X[I] += DT * vxh[i];
            Y[I] += DT * vyh[i];
            Z[I] += DT * vzh[i];
        }

        /* 3) share new positions */
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, X, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, Y, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, Z, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);

        /* 4) new acceleration */
        accel(X, Y, Z, gm_all, n, my_start, my_count, axl, ayl, azl);

        /* 5) second half kick  v <- v_half + a_new * dt/2 */
        for (int i = 0; i < my_count; ++i) {
            int I = my_start + i;
            VX[I] = vxh[i] + 0.5 * DT * axl[i];
            VY[I] = vyh[i] + 0.5 * DT * ayl[i];
            VZ[I] = vzh[i] + 0.5 * DT * azl[i];
        }

        t += DT;
    }

    if (rank == 0) {
        fclose(out);
        printf("done: %d bodies, %d steps, %.1f days -> %s\n",
               n, n_steps, YEARS * 365.25, outfile);
    }

    free(X); free(Y); free(Z);
    free(VX); free(VY); free(VZ);
    free(axl); free(ayl); free(azl);
    free(vxh); free(vyh); free(vzh);
    free(cnts); free(offs);

    MPI_Finalize();
    return 0;
}