/* -*- linux-c -*- */
/* fewbody.c

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/times.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include "fewbody.h"

fb_debug = 0;

fb_ret_t fewbody(fb_input_t input, fb_units_t units, fb_hier_t *hier, double *t, gsl_rng *rng)
{
	int i, j, k, status, done=0, forceclassify=0, restart, restep;
	long clk_tck;
	double s, slast, sstop=FB_SSTOP, tout, h=FB_H, *y, *f_quad=NULL, texpand, tnew, R[3];
	double Ei, E, Lint[3], Li[3], L[3], DeltaL[3];
	double dedt_gw_old,dedt_gw_new;
	double m,rx,ry,rz,vx,vy,vz,ax,ay,az;
	double Q_xx,Q_yy,Q_zz,Q_xy,Q_yz,Q_xz;
	double r_in_M;
	double s2, s2prev=GSL_POSINF, s2prevprev=GSL_POSINF, s2minprev=GSL_POSINF, s2max=0.0, s2min;
	struct tms firsttimebuf, currtimebuf;
	clock_t firstclock, currclock;
	FILE *quad_file=NULL;
	if(input.fb_quad){
		quad_file = fopen("Q_double_dot.dat","w"); 
		fprintf(quad_file, "#1:t(CGS) #2:Q_xx #3:Q_yy #4:Q_zz #5:Q_xy #6:Q_xz #7:Q_yz\n");
	}
	fb_hier_t phier;
	fb_ret_t retval;
	fb_nonks_params_t nonks_params;
	fb_ks_params_t ks_params;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH], logentry[FB_MAX_LOGENTRY_LENGTH];
	const gsl_odeiv_step_type *ode_type=gsl_odeiv_step_rk8pd;
	gsl_odeiv_step *ode_step;
	gsl_odeiv_control *ode_control;
	gsl_odeiv_evolve *ode_evolve;
	gsl_odeiv_system ode_sys;

	/* initialize a few things */
	fb_init_hier(hier);
	retval.iclassify = 0;
	retval.Rmin = FB_RMIN;
	retval.Rmin_i = -1;
	retval.Rmin_j = -1;
    for (k=0; k<3; k++) {
	    hier->hier[hier->hi[1]+k].Rmin = FB_RMIN;
	    hier->hier[hier->hi[1]+k].Rmin_j = -1;
    }
	retval.Nosc = 0;
	retval.DeltaE_GW = 0.;
	strncpy(logentry, input.firstlogentry, FB_MAX_LOGENTRY_LENGTH);

	/* set up the perturbation tree, initially flat */
	phier.nstarinit = hier->nstar;
	phier.nstar = hier->nstar;
	fb_malloc_hier(&phier);
	fb_init_hier(&phier);
	for (i=0; i<phier.nstar; i++) {
		fb_objcpy(&(phier.hier[phier.hi[1]+i]), &(hier->hier[hier->hi[1]+i]));
	}

	/* initialize GSL integration routine */
	if (input.ks) {
		ode_step = gsl_odeiv_step_alloc(ode_type, 8 * (hier->nstar * (hier->nstar - 1) / 2) + 1);
		ode_control = gsl_odeiv_control_y_new(input.absacc, input.relacc);
		ode_evolve = gsl_odeiv_evolve_alloc(8 * (hier->nstar * (hier->nstar - 1) / 2) + 1);
		ode_sys.function = fb_ks_func;
		ode_sys.jacobian = NULL;
		ode_sys.dimension = 8 * (hier->nstar * (hier->nstar - 1) / 2) + 1;
		ode_sys.params = &ks_params;
	} else {
		ode_step = gsl_odeiv_step_alloc(ode_type, 6 * hier->nstar);
		ode_control = gsl_odeiv_control_y_new(input.absacc, input.relacc);
		ode_evolve = gsl_odeiv_evolve_alloc(6 * hier->nstar);
		ode_sys.function = fb_nonks_func;
		ode_sys.jacobian = fb_nonks_jac;
		ode_sys.dimension = 6 * hier->nstar;
		ode_sys.params = &nonks_params;
	}

	/* set parameters for integrator */
	if (input.ks) {
		ks_params.nstar = hier->nstar;
		ks_params.kstar = ks_params.nstar*(ks_params.nstar-1)/2;
		fb_malloc_ks_params(&ks_params);
		fb_init_ks_params(&ks_params, *hier);
	} else {
		nonks_params.nstar = hier->nstar;
		fb_malloc_nonks_params(&nonks_params);
		fb_init_nonks_params(&nonks_params, *hier);
		nonks_params.PN1 = input.PN1;
		nonks_params.PN2 = input.PN2;
		nonks_params.PN25 = input.PN25;
		nonks_params.PN3 = input.PN3;
		nonks_params.PN35 = input.PN35;
		nonks_params.units = units;
	}
    if (input.PN1 == 1 || input.PN2 == 1 || input.PN25 == 1 || input.PN3 == 1 || input.PN35 == 1){
        retval.PN_ON = 1;
    }
	
	/* set the initial conditions in y_i */
	if (input.ks) {
		y = fb_malloc_vector(8*ks_params.kstar+1);
		if(input.fb_quad)
			f_quad = fb_malloc_vector(8*ks_params.kstar+1);
		y[0] = *t;
		fb_euclidean_to_ks(phier.obj, y, ks_params.nstar, ks_params.kstar);
		s = 0.0;
	} else {
		y = fb_malloc_vector(6*nonks_params.nstar);
		if(input.fb_quad)
			f_quad = fb_malloc_vector(6*nonks_params.nstar);
		fb_euclidean_to_nonks(phier.obj, y, nonks_params.nstar);
		s = *t;
	}

	/* store the initial energy and angular momentum */
	Ei = fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar) + fb_ketot(&(hier->hier[hier->hi[1]]), hier->nstar) + 
		fb_einttot(&(hier->hier[hier->hi[1]]), hier->nstar);
	dedt_gw_new = fb_dedt_gw(&(hier->hier[hier->hi[1]]), hier->nstar, units, nonks_params);
	fb_angmom(&(hier->hier[hier->hi[1]]), hier->nstar, Li);
	fb_angmomint(&(hier->hier[hier->hi[1]]), hier->nstar, Lint);
	for (i=0; i<3; i++) {
		Li[i] += Lint[i];
	}


	/* integrate along */
	fb_dprintf("fewbody: integrating...\n");
	
	done = 0;
	retval.count = 0;
	tout = *t;
	texpand = 0.0;
	clk_tck = sysconf(_SC_CLK_TCK);
	firstclock = times(&firsttimebuf);
	retval.tcpu = 0.0;
	// HVSS: things for printing trajectories
	FILE *fp;
	fp = fopen("trajectory0.dat", "w");
	fprintf(fp, "t,x0,y0,z0,vx0,vy0,vz0,x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2\n");
	while (*t < input.tstop && retval.tcpu < input.tcpustop && !done) {

		// HVSS: print trajectory information
		fprintf(fp, "%e,", (*t)*units.t);
		fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+0].x[0], hier->hier[hier->hi[1]+0].x[1], hier->hier[hier->hi[1]+0].x[2]);
		fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+0].v[0], hier->hier[hier->hi[1]+0].v[1], hier->hier[hier->hi[1]+0].v[2]);
		fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+1].x[0], hier->hier[hier->hi[1]+1].x[1], hier->hier[hier->hi[1]+1].x[2]);
		fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+1].v[0], hier->hier[hier->hi[1]+1].v[1], hier->hier[hier->hi[1]+1].v[2]);
		fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+2].x[0], hier->hier[hier->hi[1]+2].x[1], hier->hier[hier->hi[1]+2].x[2]);
		fprintf(fp, "%e,%e,%e\n", hier->hier[hier->hi[1]+2].v[0], hier->hier[hier->hi[1]+2].v[1], hier->hier[hier->hi[1]+2].v[2]);
		//fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", (*t)*units.t, hier->obj[0]->x[0], hier->obj[0]->x[1], hier->obj[0]->x[2], hier->obj[0]->v[0], hier->obj[0]->v[1], hier->obj[0]->v[2], hier->obj[1]->x[0], hier->obj[1]->x[1], hier->obj[1]->x[2], hier->obj[1]->v[0], hier->obj[1]->v[1], hier->obj[1]->v[2], hier->obj[2]->x[0], hier->obj[2]->x[1], hier->obj[2]->x[2], hier->obj[2]->v[0], hier->obj[2]->v[1], hier->obj[2]->v[2]);

		if (hier->nstarinit == 2 && hier->nstar == 2 && hier->nobj == 1) {
			fb_upsync(&(hier->hier[hier->hi[2]+0]), *t);
		}
		
		/* take one step */
		slast = s;
		status = gsl_odeiv_evolve_apply(ode_evolve, ode_control, ode_step, &ode_sys, &s, sstop, &h, y);
		if (status != GSL_SUCCESS) {
			break;
		}

		/* set objects' positions and velocities in phier */
		if (input.ks) {
			tnew = y[0];
			fb_ks_to_euclidean(y, phier.obj, ks_params.nstar, ks_params.kstar);
		} else {
			tnew = s;
			fb_nonks_to_euclidean(y, phier.obj, nonks_params.nstar);
		}

		/* restart: do we need to restart the integrator?
		   restep: do we need to repeat the last integration step?
		   forceclassify: do we need to force a call to the classify() routine? */
		restart = 0;
		restep = 0;
		forceclassify = 0;

		/* see if we need to expand or collapse the perturbation hierarchy */
		if (fb_expand(&phier, tnew, input.tidaltol)) {
			texpand = tnew;
			s = slast;
			restart = 1;
			restep = 1;
			for (i=0; i<hier->nstar; i++) {
				for (k=0; k<3; k++) {
					phier.hier[phier.hi[1]+i].x[k] = hier->hier[hier->hi[1]+i].x[k];
					phier.hier[phier.hi[1]+i].v[k] = hier->hier[hier->hi[1]+i].v[k];
				}
			}
			fb_elkcirt(&phier, *t);
		} else if (tnew >= texpand && fb_collapse(&phier, tnew, input.tidaltol, input.speedtol, units, nonks_params, input)) {
			*t = tnew;
			/* if there is only one object, then it's stable---force classify() */
			if (phier.nobj == 1) {
				forceclassify = 1;
			} else {
				restart = 1;
			}
		} else {
			*t = tnew;
		}
		
		/* if we're not repeating the previous integration step, then do physics */
		if (!restep) {
			/* trickle down so updated information is in hier */
			fb_trickle(&phier, *t);
			for (i=0; i<hier->nstar; i++) {
				for (k=0; k<3; k++) {
					hier->hier[hier->hi[1]+i].x[k] = phier.hier[phier.hi[1]+i].x[k];
					hier->hier[hier->hi[1]+i].v[k] = phier.hier[phier.hi[1]+i].v[k];
				}
			}

			/* Calculate the energy lost due to gravitational-wave emission.
			 * TODO: Just a trapezoid rule integration, but for any inspirals
			 * dedt will be strictly negative and decreasing, meaning the errors
			 * may be significant and systematic.  Think about using a smarter integrator */
			    if(retval.PN_ON){
				dedt_gw_old = dedt_gw_new;
				dedt_gw_new = fb_dedt_gw(&(hier->hier[hier->hi[1]]), hier->nstar, units, nonks_params);
				retval.DeltaE_GW += (dedt_gw_old + dedt_gw_new) * (s - slast) / 2.;
			    }

			/* update Rmin (closest approach) */
			for (i=0; i<hier->nstar-1; i++) {
				for (j=i+1; j<hier->nstar; j++) {
					for (k=0; k<3; k++) {
						R[k] = hier->hier[hier->hi[1]+i].x[k] - hier->hier[hier->hi[1]+j].x[k];
					}
	
                    /* For PN terms, compute a,e for any possible pairs that may soon merge */
                    if(retval.PN_ON) {
                        r_in_M = fb_compute_distance_in_M(&hier->hier[hier->hi[1]+i],&hier->hier[hier->hi[1]+j], units);
        
                        if (r_in_M < 510) {
                            fb_check_ecc_for_inspiral(&hier->hier[hier->hi[1]+i],&hier->hier[hier->hi[1]+j], r_in_M, units);
                        }
                    }
	
					if (fb_mod(R) < retval.Rmin) {
						retval.Rmin = fb_mod(R);
						retval.Rmin_i = i;
						retval.Rmin_j = j;
                    }

                    /* HVSS: checks to see if either of the star-specific Rmins need to be updated */
					if (fb_mod(R) < hier->hier[hier->hi[1]+i].Rmin) {
						hier->hier[hier->hi[1]+i].Rmin = fb_mod(R);
						hier->hier[hier->hi[1]+i].Rmin_j = j;
                    }
					if (fb_mod(R) < hier->hier[hier->hi[1]+j].Rmin) {
						hier->hier[hier->hi[1]+j].Rmin = fb_mod(R);
						hier->hier[hier->hi[1]+j].Rmin_j = i;
                    }
				}
			}
			
			/* do physical collisions */
			if (fb_collide(hier, input.fexp, units, rng, input.BH_REFF, input)) {
				/* initialize phier to a flat tree */
				phier.nstar = hier->nstar;
				fb_init_hier(&phier);
				for (i=0; i<phier.nstar; i++) {
					fb_objcpy(&(phier.hier[phier.hi[1]+i]), &(hier->hier[hier->hi[1]+i]));
				}
				
				/* if there is only one object, then it's stable---force classify() */
				if (phier.nobj == 1) {
					restart = 0;
					forceclassify = 1;
				} else {
					restart = 1;
				}
			}
			
			/* check for oscillations in s^2 */
			/* first calculate s^2 */
			s2 = 0.0;
			for (i=0; i<hier->nstar-1; i++) {
				for (j=i+1; j<hier->nstar; j++) {
					for (k=0; k<3; k++) {
						s2 += fb_sqr(hier->hier[hier->hi[1]+i].x[k] - hier->hier[hier->hi[1]+j].x[k]);
					}
				}
			}
			/* local min */
			if (s2 > s2prev && s2prev < s2prevprev) {
				s2min = s2prev;
				/* increment Nosc if there is a valid oscillation */
				if (s2max >= 2.0*s2minprev && s2max >= 2.0*s2min) {
					retval.Nosc++;
					fb_dprintf("fewbody: oscillation in s^2 detected: s2max/s2minprev=%g s2max/s2min=%g: Nosc=%d\n", s2max/s2minprev, s2max/s2min, retval.Nosc);
					s2minprev = s2min;
				} else {
					fb_dprintf("fewbody: oscillation in s^2 not big enough: s2max/s2minprev=%g s2max/s2min=%g\n", s2max/s2minprev, s2max/s2min);
					/* that pesky first min... */
					if (s2max == 0.0) {
						s2minprev = s2min;
					}
				}
			}
			/* local max: just set value and wait for next min */
			if (s2 < s2prev && s2prev > s2prevprev) {
				s2max = s2prev;
			}
			/* set previous values */
			s2prevprev = s2prev;
			s2prev = s2;

			/* see if we're done */
			if (retval.count % input.ncount == 0 || forceclassify) {
				status = fb_classify(hier, *t, input.tidaltol, input.speedtol, units, input);
				retval.iclassify++;
				fb_dprintf("fewbody: current status:  t=%.6g  %s  (%s)\n",
					   *t, fb_sprint_hier(*hier, string1),
					   fb_sprint_hier_hr(*hier, string2));
				snprintf(&(logentry[strlen(logentry)]), FB_MAX_LOGENTRY_LENGTH-strlen(logentry),
					 "  current status:  t=%.6g  %s  (%s)\n", *t, fb_sprint_hier(*hier, string1),
					 fb_sprint_hier_hr(*hier, string2));
				if (status) {
					done = 1;
				}
			}
			
			/* print stuff if necessary */
			if (input.Dflag == 1 && (*t >= tout || done)) {
				tout = *t + input.dt;
				fb_print_story(&(hier->hier[hier->hi[1]]), hier->nstar, *t, logentry);
			} else if (input.Dflag == 2){
				fb_print_story(&(hier->hier[hier->hi[1]]), hier->nstar, *t, logentry);
			}
		}
		
		/* restart integrator if necessary */
		if (restart) {
			fb_dprintf("fewbody: restarting integrator: nobj=%d count=%ld\n", phier.nobj, retval.count);
			fb_free_vector(y);
			if (input.ks) {
				fb_free_ks_params(ks_params);
				ks_params.nstar = phier.nobj;
				ks_params.kstar = ks_params.nstar*(ks_params.nstar-1)/2;
				fb_malloc_ks_params(&ks_params);
				fb_init_ks_params(&ks_params, phier);
				
				y = fb_malloc_vector(8*ks_params.kstar+1);
				if(input.fb_quad)
					f_quad = fb_malloc_vector(8*ks_params.kstar+1);
				y[0] = *t;
				fb_euclidean_to_ks(phier.obj, y, ks_params.nstar, ks_params.kstar);
			} else {
				fb_free_nonks_params(nonks_params);
				nonks_params.nstar = phier.nobj;
				fb_malloc_nonks_params(&nonks_params);
				fb_init_nonks_params(&nonks_params, phier);
				nonks_params.PN1 = input.PN1;
				nonks_params.PN2 = input.PN2;
				nonks_params.PN25 = input.PN25;
				nonks_params.PN3 = input.PN3;
				nonks_params.PN35 = input.PN35;
				nonks_params.units = units;
				
				y = fb_malloc_vector(6*nonks_params.nstar);
				if(input.fb_quad)
					f_quad = fb_malloc_vector(6*nonks_params.nstar);
				fb_euclidean_to_nonks(phier.obj, y, nonks_params.nstar);
			}
			
			/* and re-allocate integrator */
			gsl_odeiv_evolve_free(ode_evolve);
			gsl_odeiv_control_free(ode_control);
			gsl_odeiv_step_free(ode_step);
			
			/* re-initialize integrator */
			if (input.ks) {
				ode_step = gsl_odeiv_step_alloc(ode_type, 8*ks_params.kstar+1);
				ode_control = gsl_odeiv_control_y_new(input.absacc, input.relacc);
				ode_evolve = gsl_odeiv_evolve_alloc(8*ks_params.kstar+1);
				ode_sys.dimension = 8*ks_params.kstar+1;
			} else {
				ode_step = gsl_odeiv_step_alloc(ode_type, 6*nonks_params.nstar);
				ode_control = gsl_odeiv_control_y_new(input.absacc, input.relacc);
				ode_evolve = gsl_odeiv_evolve_alloc(6*nonks_params.nstar);
				ode_sys.dimension = 6*nonks_params.nstar;
			}
		}

		if(!restep && !restart && input.fb_quad){
			fb_nonks_func(*t, y, f_quad, &nonks_params);
			Q_xx = 0;
			Q_yy = 0;
			Q_zz = 0;
			Q_xy = 0;
			Q_yz = 0;
			Q_xz = 0; 

			for (i=0; i<nonks_params.nstar; i++){ 
				m = nonks_params.m[i] * units.m;
				rx = y[i*6 + 0] * units.l;
				ry = y[i*6 + 1] * units.l;
				rz = y[i*6 + 2] * units.l;
				vx = y[i*6 + 3] * units.l / units.t;
				vy = y[i*6 + 4] * units.l / units.t;
				vz = y[i*6 + 5] * units.l / units.t;
				ax = f_quad[i*6 + 3] * units.l / units.t / units.t / m;
				ay = f_quad[i*6 + 4] * units.l / units.t / units.t / m;
				az = f_quad[i*6 + 5] * units.l / units.t / units.t / m;

				Q_xx += m*(vx*vx + rx*ax  - 2*vy*vy - 2*ry*ay - 2*vz*vz - 2*rz*az)/3.;
				Q_yy += m*(vy*vy + ry*ay  - 2*vx*vx - 2*rx*ax - 2*vz*vz - 2*rz*az)/3.;
				Q_zz += m*(vz*vz + rz*az  - 2*vy*vy - 2*ry*ay - 2*vx*vx - 2*rx*ax)/3.;
				Q_xy += m*(ax*ry + 2*vx*vy + rx*ay);
				Q_yz += m*(ay*rz + 2*vy*vz + ry*az);
				Q_xz += m*(ax*rz + 2*vx*vz + rx*az);
			}

			fprintf(quad_file, "%.18e %.8e %.8e %.8e %.8e %.8e %.8e \n",(*t)*units.t,Q_xx,Q_yy,Q_zz,Q_xy,Q_yz,Q_xz);
		}

			

		/* update variables that change on every integration step */
		retval.count++;
		currclock = times(&currtimebuf);
		retval.tcpu = ((double) (currtimebuf.tms_utime + currtimebuf.tms_stime - firsttimebuf.tms_utime - firsttimebuf.tms_stime))/((double) clk_tck);
	}
	// HVSS: write last step, and close file
	fprintf(fp, "%e,", (*t)*units.t);
	fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+0].x[0], hier->hier[hier->hi[1]+0].x[1], hier->hier[hier->hi[1]+0].x[2]);
	fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+0].v[0], hier->hier[hier->hi[1]+0].v[1], hier->hier[hier->hi[1]+0].v[2]);
	fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+1].x[0], hier->hier[hier->hi[1]+1].x[1], hier->hier[hier->hi[1]+1].x[2]);
	fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+1].v[0], hier->hier[hier->hi[1]+1].v[1], hier->hier[hier->hi[1]+1].v[2]);
	fprintf(fp, "%e,%e,%e,", hier->hier[hier->hi[1]+2].x[0], hier->hier[hier->hi[1]+2].x[1], hier->hier[hier->hi[1]+2].x[2]);
	fprintf(fp, "%e,%e,%e\n", hier->hier[hier->hi[1]+2].v[0], hier->hier[hier->hi[1]+2].v[1], hier->hier[hier->hi[1]+2].v[2]);
	//fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", (*t)*units.t, hier->obj[0]->x[0], hier->obj[0]->x[1], hier->obj[0]->x[2], hier->obj[0]->v[0], hier->obj[0]->v[1], hier->obj[0]->v[2], hier->obj[1]->x[0], hier->obj[1]->x[1], hier->obj[1]->x[2], hier->obj[1]->v[0], hier->obj[1]->v[1], hier->obj[1]->v[2], hier->obj[2]->x[0], hier->obj[2]->x[1], hier->obj[2]->x[2], hier->obj[2]->v[0], hier->obj[2]->v[1], hier->obj[2]->v[2]);
	fclose(fp);
	
	/* do final classification */
	retval.retval = fb_classify(hier, *t, input.tidaltol, input.speedtol, units, input);
	retval.iclassify++;
	fb_dprintf("fewbody: current status:  t=%.6g  %s  (%s)\n",
		   *t, fb_sprint_hier(*hier, string1),
		   fb_sprint_hier_hr(*hier, string2));
	snprintf(&(logentry[strlen(logentry)]), FB_MAX_LOGENTRY_LENGTH-strlen(logentry),
		 "  current status:  t=%.6g  %s  (%s)\n", *t, fb_sprint_hier(*hier, string1),
		 fb_sprint_hier_hr(*hier, string2));
	
	/* print final story */
	if (input.Dflag == 1) {
		fb_print_story(&(hier->hier[hier->hi[1]]), hier->nstar, *t, logentry);
	}
	
	fb_dprintf("fewbody: final: phier.nobj = %d\n", phier.nobj);

	E = fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar) + fb_ketot(&(hier->hier[hier->hi[1]]), hier->nstar) + 
		fb_einttot(&(hier->hier[hier->hi[1]]), hier->nstar);
	fb_angmom(&(hier->hier[hier->hi[1]]), hier->nstar, L);
	fb_angmomint(&(hier->hier[hier->hi[1]]), hier->nstar, Lint);
	for (i=0; i<3; i++) {
		L[i] += Lint[i];
		DeltaL[i] = L[i] - Li[i];
	}

	/* free GSL stuff */
	gsl_odeiv_evolve_free(ode_evolve);
	gsl_odeiv_control_free(ode_control);
	gsl_odeiv_step_free(ode_step);

	/* free our own stuff */
	fb_free_vector(y);
	fb_free_hier(phier);

	if (input.ks) {
		fb_free_ks_params(ks_params);
	} else {
		fb_free_nonks_params(nonks_params);
	}

	/* done! */
	retval.DeltaE = E-Ei;
	retval.DeltaEfrac = E/Ei-1.0;
	retval.DeltaL = fb_mod(DeltaL);
	retval.DeltaLfrac = fb_mod(DeltaL)/fb_mod(Li);
   	if (retval.PN_ON)
		retval.DeltaE_GWfrac = retval.DeltaE_GW/Ei; /* technically this should be the pN Ei, but any scattering should start off in the Newtonian regime anyway...*/
	if(input.fb_quad)
		fclose(quad_file);
	return(retval);
}
