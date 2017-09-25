#!/usr/bin/env python
#
# EOS_Tables: manages EOS tables for THC
# Loosely based on Filippo Galeazzi's Matlab scripts
# Copyright (C) 2016, David Radice <dradice@caltech.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import h5py
import numpy as np
import sys
import os

sys.path.append(os.path.dirname(sys.argv[0]) + "/../modules")
import betaeq
from eos_bar_table import COLWIDTH, EOS_Table
import unitconv as ut

INFO = "INFO".ljust(COLWIDTH)

# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser()
# -----------------------------------------------------------------------------
parser.add_argument("-y", "--hydro", dest="hydro", required=True,
    help="Hydro part of the EOS table")
parser.add_argument("-w", "--weak", dest="weak", required=True,
    help="Weak part of the EOS table")
parser.add_argument("-o", "--output", dest="output", required=True,
    help="Output base file name (will be overwritten)")
parser.add_argument("-t", "--temperature", dest="temp", type=float,
    help="Slice temperature in MeV (defaults to the minimum on the table)")
parser.add_argument("--rm-radiation", dest="rm_radiation", action="store_true",
    help="Remove the radiation pressure part of the EOS")
parser.add_argument("--attach-poly", dest="attach_poly", action="store_true",
    help="Attach a polytrope at low density")
parser.add_argument("--resample", dest="resample", type=int, default=-1,
    help="Resample the EOS table, if negative don't resample (default: -1)")
parser.add_argument("-l", "--lorene", dest="lorene", action="store_true",
    help="Output the EOS slice in LORENE format")
parser.add_argument("-p", "--pizza", dest="pizza", action="store_true",
    help="Output the EOS slice in PIZZA format")
parser.add_argument("-r", "--rns", dest="rns", action="store_true",
    help="Output the EOS slice in RNS format")
args = parser.parse_args()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Read the table into memory, compute quantities of interest
# -----------------------------------------------------------------------------
table = {}
print(INFO + "hydro table: {}".format(args.hydro))
dfile = h5py.File(args.hydro, "r")
table["rho"] = np.array(dfile["density"])
table["temp"] = np.array(dfile["temperature"])
table["ye"]  = np.array(dfile["ye"])
table["cs2"] = np.array(dfile["cs2"])
table["eps"] = np.array(dfile["internalEnergy"])
table["press"] = np.array(dfile["pressure"])
table["mass_factor"] = np.array(dfile["mass_factor"])
del dfile

print(INFO + "weak table: {}".format(args.weak))
dfile = h5py.File(args.weak, "r")
table["mu_e"] = np.array(dfile["mu_e"])
table["mu_p"] = np.array(dfile["mu_p"])
table["mu_n"] = np.array(dfile["mu_n"])
del dfile

mass_factor_cgs = table["mass_factor"]*(ut.MEV_CGS/(ut.C_CGS**2))
nb = table["rho"]/(mass_factor_cgs/(ut.FM_CGS**3))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Compute beta equilibrium slice
# -----------------------------------------------------------------------------
if args.temp is not None:
    itemp = np.argmin(np.abs(args.temp - table["temp"]))
else:
    itemp = 0
temp = table["temp"][itemp]
print(INFO + "making beta-equilibrium slice at T = {} MeV".format(temp))

ye_beta = np.empty_like(table["rho"])
rho0_beta = table["rho"]
ye_beta[:] = np.NAN
for inb in range(nb.shape[0]):
    ye_beta[inb] = betaeq.find_beta_eq(table["ye"],
            table["mu_e"][:,itemp,inb], table["mu_n"][:,itemp,inb],
            table["mu_p"][:,itemp,inb])
assert np.all(np.isfinite(ye_beta))

cs2_beta = np.empty_like(ye_beta)
eps_beta = np.empty_like(ye_beta)
press_beta = np.empty_like(ye_beta)
for inb in range(nb.shape[0]):
    cs2_beta[inb] = betaeq.interp_f_of_ye(table["ye"],
            table["cs2"][:,itemp,inb])(ye_beta[inb])
    eps_beta[inb] = betaeq.interp_f_of_ye(table["ye"],
            table["eps"][:,itemp,inb])(ye_beta[inb])
    press_beta[inb] = betaeq.interp_f_of_ye(table["ye"],
            table["press"][:,itemp,inb])(ye_beta[inb])
assert np.all(np.isfinite(cs2_beta))
assert np.all(np.isfinite(eps_beta))
assert np.all(np.isfinite(press_beta))

# Clean soundspeed
cs2_beta[cs2_beta < 0.0] = 0.0
cs2_beta[cs2_beta > 1.0] = 1.0

# Total energy density in g/cc
rho_beta = nb * (mass_factor_cgs/(ut.FM_CGS**3)) * \
    (1.0 + (eps_beta/(ut.C_CGS**2)))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Remove the radiation pressure / energy
# -----------------------------------------------------------------------------
if args.rm_radiation:
    print(INFO + "removing radiation pressure")
    rad_press = 1.0/3.0 * ut.RAD_CGS * (temp * ut.MEV_CGS/ut.KB_CGS)**4
    rad_rho = rad_press * 3.0 / (ut.C_CGS**2)

    press_beta -= rad_press
    rho_beta -= rad_rho

    eps_beta = (rho_beta - rho0_beta)/rho0_beta * ut.C_CGS**2

    assert np.all(press_beta > 0)
    assert np.all(rho_beta > 0)
    assert np.all(eps_beta > 0)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Create barotropic EOS table object
# -----------------------------------------------------------------------------
mb_SI = (table["mass_factor"]*ut.MEV_SI)/(ut.C_SI**2)
mb_PU = mb_SI / ut.PIZZA_UNITS.mass
c_PU  = ut.C_SI / ut.PIZZA_UNITS.velocity
uc    = ut.CGS_UNITS / ut.PIZZA_UNITS
eos_slice = EOS_Table(
    rho0_beta * uc.density,
    eps_beta/(ut.C_CGS**2),
    press_beta * uc.pressure,
    efr  = ye_beta,
    temp = np.ones_like(ye_beta) * temp,
    mbar = mb_PU,
    csnd = np.sqrt(cs2_beta)*c_PU,
    name = args.output)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Attach polytrope at low density and resample
# -----------------------------------------------------------------------------
if args.attach_poly:
    print(INFO + "attaching polytrope")
    eos_slice = eos_slice.make_restmass_natural(3)
    eos_slice = eos_slice.make_poly_compatible(3)
    if args.resample > 0:
        npoints = args.resample
    else:
        npoints = eos_slice.rmd.shape[0] + 50
    eos_slice = eos_slice.attach_poly(1*uc.density, npoints)
    eos_slice = eos_slice.make_adiabatic()
elif args.resample > 0:
    print(INFO + "resampling the EOS slice")
    eos_slice = eos_slice.resample_geom(args.resample)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Output the data in ASCII format
# -----------------------------------------------------------------------------
print(INFO + "exporting tables")
if args.lorene:
    eos_slice.save_lorene("{}.lorene".format(args.output))
if args.pizza:
    eos_slice.save_pizza("{}.pizza".format(args.output))
if args.rns:
    eos_slice.save_rns("{}.rns".format(args.output))
# -----------------------------------------------------------------------------
