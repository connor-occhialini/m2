#include <stdio.h>
#include <math.h>

#include "ising_datastructure.h"
#include "ising.h"
#include "metropolis.h"
#include "progressbar.h"
#include "ising-demo-metropolis.h"
#include "plotfun.h"

int main (int argc, char **argv)
{
    int i;

    struct ising model;

    double tmin=3.5, tmax=3.8, dt,t=tmin;       //temperature range           
    double h = 0.;                  	        //magnetic field (held at zero)
    double j = 1.;                 		//coupling strength 
    int nt = 8;					//total number of temperatures to compute
    dt = (tmax-tmin)/(nt-1);			//temperature increment
    int lx = 32;
    int ly = lx;
    double specheat[nt];			//specific heat value array
    double temperat[nt];			//temperature value array
    
    for(int v = 0; v < nt; v++)			//temperature array initialization
    {
        temperat[v] = tmin + v*dt;
    }


    int mcsteps = 1024 * 128;        		//total number of steps

    if( ising_init (&model, lx, ly, j, h, t) != 0)
    {
        fprintf(stderr, "%s - cannot allocate spin array\n", argv[0]);
        argc = 1;
        exit(1);
    }

    print_header (&model, stderr);

    int therm_steps = mcsteps / 4;
    int report = therm_steps / 32;
    
    for(i = 0; i < nt; i++)				//loop over temperatures
    {
    for (int s = 1; s <= therm_steps; s++)
    {
        (void) one_metropolis_step_per_spin (&model);
        if (s % report == 0)
        {
            text_progress_bar (s, therm_steps, 32, "Thermalization", stderr);
        }
    }
    fprintf (stderr, " Done\n");

    report = mcsteps / 32;
    double mav = 0., m2av = 0., eav = 0., e2av = 0.;
    double cv, xi;

    for (int s = 1; s <= mcsteps; s++)
    {
        (void) one_metropolis_step_per_spin (&model);
        double m = magnetization_per_spin (&model);
        double e = energy_per_spin (&model);

        mav += m;
        m2av += m * m;
        eav += e;
        e2av += e * e;

        if (s % report == 0)
        {
            text_progress_bar (s, mcsteps, 32, "Production", stderr);
        }

    }
    fprintf (stderr, " Done\n\n");

    mav /= mcsteps;
    eav /= mcsteps;
    m2av /= mcsteps;
    e2av /= mcsteps;

    // specific heat dE/dt
    cv = (e2av - eav*eav)/(t*t);
    specheat[i] = cv;

    // magnetic susceptibility
    xi = (m2av - mav*mav)/t;

    printf ("%.3f   %.3f     % f  %f     % f  %f\n", t, h, mav, 
        xi, eav, cv);
    t += dt;
    ising_reinit (&model,h,t);
    }
    ising_free (model);
    pfun(temperat, specheat, nt);			//function plots specific heat vs. temperature
    return 0;
}

void print_header (struct ising *model, FILE *fp)
{
    fprintf (fp, "\n");
    fprintf (fp, " Two-dimensional Ising model - Metropolis simulation\n");
    fprintf (fp, " ---------------------------------------------------\n");
    fprintf (fp, " T/J = %.3f  H/J = %.3f  %dx%d\n\n", model->T/model->J, 
        model->H/model->J, model->Lx, model->Ly);
}

