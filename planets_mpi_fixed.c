/*
 * planets_mpi.c
 * Parallel N-body simulation using MPI and velocity Verlet.
 *
 * Input formats supported:
 *   10 columns: name,parent,gm_km3_s2,gm_au3_d2,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud
 *    9 columns: name,parent,gm_km3_s2,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud
 *
 * Build:
 *   mpicc -O2 -Wall -Wextra -o planets_mpi planets_mpi.c -lm
 *
 * Run:
 *   mpirun -np 4 ./planets_mpi initial.csv traj.csv
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define MAX_BODIES 128
#define MAX_NAME   64
#define AU_KM      149597870.7
#define DAY_S      86400.0
#define YEARS      5.0
#define DT         0.0001       /* days per integration step */
#define SAVE_EVERY 10000        /* 0.0001 * 10000 = 1 day */

typedef struct {
    char name[MAX_NAME];
    double gm;                  /* AU^3 / day^2 */
    double x, y, z;             /* AU */
    double vx, vy, vz;          /* AU / day */
} Body;

static Body all[MAX_BODIES];
static double gm_all[MAX_BODIES];
static int n = 0;

/* Read the initial conditions on rank 0. */
static int read_initial(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error: cannot open input file '%s'\n", filename);
        return -1;
    }

    char line[1024];
    int skipped_header = 0;

    while (fgets(line, sizeof line, f) != NULL) {
        /* Skip blank lines and comments. */
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;

        /* Skip the CSV header (the first non-comment line). */
        if (!skipped_header) {
            skipped_header = 1;
            continue;
        }

        int commas = 0;
        for (const char *p = line; *p != '\0'; ++p) {
            if (*p == ',')
                ++commas;
        }

        char name[MAX_NAME] = {0};
        char parent[MAX_NAME] = {0};
        double gm_km = 0.0;
        double gm_au_unused = 0.0;
        double x, y, z, vx, vy, vz;
        int fields_read;

        if (commas >= 9) {
            /* 10-column format: both km^3/s^2 and AU^3/day^2 GM values. */
            fields_read = sscanf(
                line,
                "%63[^,],%63[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                name, parent, &gm_km, &gm_au_unused,
                &x, &y, &z, &vx, &vy, &vz
            );

            if (fields_read != 10)
                continue;
        } else {
            /* 9-column format: GM in km^3/s^2 only. */
            fields_read = sscanf(
                line,
                "%63[^,],%63[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                name, parent, &gm_km, &x, &y, &z, &vx, &vy, &vz
            );

            if (fields_read != 9)
                continue;
        }

        if (n >= MAX_BODIES) {
            fprintf(stderr, "Warning: maximum of %d bodies reached; extra rows ignored.\n",
                    MAX_BODIES);
            break;
        }

        Body *b = &all[n];
        snprintf(b->name, sizeof b->name, "%s", name);

        /* Convert GM from km^3/s^2 to AU^3/day^2. */
        b->gm = gm_km * DAY_S * DAY_S /
                (AU_KM * AU_KM * AU_KM);

        b->x = x;
        b->y = y;
        b->z = z;
        b->vx = vx;
        b->vy = vy;
        b->vz = vz;

        ++n;
    }

    fclose(f);
    return n;
}

/* Calculate acceleration for the subset of bodies owned by this rank. */
static void accel(const double *X, const double *Y, const double *Z,
                  const double *GM, int ntot,
                  int start, int count,
                  double *AX, double *AY, double *AZ)
{
    for (int i = 0; i < count; ++i) {
        const int I = start + i;
        double ax = 0.0, ay = 0.0, az = 0.0;

        for (int J = 0; J < ntot; ++J) {
            if (J == I)
                continue;

            const double dx = X[J] - X[I];
            const double dy = Y[J] - Y[I];
            const double dz = Z[J] - Z[I];

            /*
             * Small softening term reduces singular behavior at extremely
             * small separations; it is not a physical collision model.
             */
            const double r2 = dx * dx + dy * dy + dz * dz + 1e-12;
            const double k = GM[J] / (r2 * sqrt(r2));

            ax += k * dx;
            ay += k * dy;
            az += k * dz;
        }

        AX[i] = ax;
        AY[i] = ay;
        AZ[i] = az;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const char *infile = (argc > 1) ? argv[1] : "initial.csv";
    const char *outfile = (argc > 2) ? argv[2] : "traj.csv";

    /* Only rank 0 reads the CSV file. */
    if (rank == 0) {
        if (read_initial(infile) < 0)
            MPI_Abort(MPI_COMM_WORLD, 1);

        if (n <= 0) {
            fprintf(stderr, "Error: no valid bodies found in '%s'\n", infile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* First broadcast the number of bodies. */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*
     * IMPORTANT: each MPI process has its own memory.
     * Broadcast the complete Body array, not just n.
     */
    MPI_Bcast(all, (int)(n * sizeof(Body)), MPI_BYTE, 0, MPI_COMM_WORLD);

    /* Distribute body indices as evenly as possible among ranks. */
    const int base = n / size;
    const int extra = n % size;
    const int my_start = rank * base + (rank < extra ? rank : extra);
    const int my_count = base + (rank < extra ? 1 : 0);

    for (int i = 0; i < n; ++i)
        gm_all[i] = all[i].gm;

    /* Global state: every rank stores a full copy of positions and velocities. */
    double *X  = calloc((size_t)n, sizeof *X);
    double *Y  = calloc((size_t)n, sizeof *Y);
    double *Z  = calloc((size_t)n, sizeof *Z);
    double *VX = calloc((size_t)n, sizeof *VX);
    double *VY = calloc((size_t)n, sizeof *VY);
    double *VZ = calloc((size_t)n, sizeof *VZ);

    /* Allocate at least one element even if this rank owns zero bodies. */
    const size_t local_size = (my_count > 0) ? (size_t)my_count : 1U;
    double *axl = calloc(local_size, sizeof *axl);
    double *ayl = calloc(local_size, sizeof *ayl);
    double *azl = calloc(local_size, sizeof *azl);
    double *vxh = calloc(local_size, sizeof *vxh);
    double *vyh = calloc(local_size, sizeof *vyh);
    double *vzh = calloc(local_size, sizeof *vzh);

    int *cnts = malloc((size_t)size * sizeof *cnts);
    int *offs = malloc((size_t)size * sizeof *offs);

    if (!X || !Y || !Z || !VX || !VY || !VZ ||
        !axl || !ayl || !azl || !vxh || !vyh || !vzh ||
        !cnts || !offs) {
        fprintf(stderr, "Rank %d: memory allocation failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int r = 0; r < size; ++r) {
        cnts[r] = base + (r < extra ? 1 : 0);
        offs[r] = r * base + (r < extra ? r : extra);
    }

    /* Initialize the full state from the broadcast Body array. */
    for (int i = 0; i < n; ++i) {
        X[i] = all[i].x;
        Y[i] = all[i].y;
        Z[i] = all[i].z;
        VX[i] = all[i].vx;
        VY[i] = all[i].vy;
        VZ[i] = all[i].vz;
    }

    /* Initial accelerations for the bodies owned by this rank. */
    accel(X, Y, Z, gm_all, n, my_start, my_count, axl, ayl, azl);

    FILE *out = NULL;
    if (rank == 0) {
        out = fopen(outfile, "w");
        if (out == NULL) {
            fprintf(stderr, "Error: cannot open output file '%s'\n", outfile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fprintf(out, "t_days,name,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud\n");
    }

    const int n_steps = (int)(YEARS * 365.25 / DT + 0.5);
    double t = 0.0;

    for (int step = 0; step <= n_steps; ++step) {
        /* Save a snapshot. All state arrays must be synchronized first. */
        if (step % SAVE_EVERY == 0) {
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                           VX, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                           VY, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                           VZ, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);

            if (rank == 0) {
                for (int i = 0; i < n; ++i) {
                    fprintf(out,
                            "%.3f,%s,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                            t, all[i].name,
                            X[i], Y[i], Z[i], VX[i], VY[i], VZ[i]);
                }
            }
        }

        if (step == n_steps)
            break;

        /* 1. Half kick: v_half = v + a * dt/2. */
        for (int i = 0; i < my_count; ++i) {
            const int I = my_start + i;
            vxh[i] = VX[I] + 0.5 * DT * axl[i];
            vyh[i] = VY[I] + 0.5 * DT * ayl[i];
            vzh[i] = VZ[I] + 0.5 * DT * azl[i];
        }

        /* 2. Drift: x = x + v_half * dt. */
        for (int i = 0; i < my_count; ++i) {
            const int I = my_start + i;
            X[I] += DT * vxh[i];
            Y[I] += DT * vyh[i];
            Z[I] += DT * vzh[i];
        }

        /* 3. Share updated positions with every process. */
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       X, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       Y, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       Z, cnts, offs, MPI_DOUBLE, MPI_COMM_WORLD);

        /* 4. Calculate new accelerations. */
        accel(X, Y, Z, gm_all, n, my_start, my_count, axl, ayl, azl);

        /* 5. Second half kick: v = v_half + a_new * dt/2. */
        for (int i = 0; i < my_count; ++i) {
            const int I = my_start + i;
            VX[I] = vxh[i] + 0.5 * DT * axl[i];
            VY[I] = vyh[i] + 0.5 * DT * ayl[i];
            VZ[I] = vzh[i] + 0.5 * DT * azl[i];
        }

        t += DT;
    }

    if (rank == 0) {
        fclose(out);
        printf("Done: %d bodies, %d steps, %.1f days -> %s\n",
               n, n_steps, YEARS * 365.25, outfile);
    }

    free(X);
    free(Y);
    free(Z);
    free(VX);
    free(VY);
    free(VZ);
    free(axl);
    free(ayl);
    free(azl);
    free(vxh);
    free(vyh);
    free(vzh);
    free(cnts);
    free(offs);

    MPI_Finalize();
    return 0;
}
