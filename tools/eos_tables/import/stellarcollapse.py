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
import datetime
import h5py
import numpy as np
import sys
import os

sys.path.append(os.path.dirname(sys.argv[0]) + "/../modules")
import unitconv as ut
import tabentries

# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser()
# -----------------------------------------------------------------------------
parser.add_argument("table", metavar="table.h5", nargs=1,
        help="Input table in stallarcollapse.org format")
parser.add_argument("-n", "--name", dest="name", required=True,
        help="Name of the EOS")
parser.add_argument("-m", "--baryon-mass", dest="mass_factor", type=float,
        help="Fiducial baryon mass in MeV. "\
            "If not specified, it is computed internally")
parser.add_argument("--rho-max", dest="rho_max", type=float, default=0.0,
        help="Exclude regions with densities larger than rho_max")
args = parser.parse_args()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Read the table into memory
# -----------------------------------------------------------------------------
dfile = h5py.File(args.table[0], "r")
table = {}
for k in dfile.keys():
    table[k] = np.array(dfile[k])
del dfile

# Remove high density region
if args.rho_max > 0:
  rho = 10.0**table["logrho"]
  irmax = np.searchsorted(rho, args.rho_max, side='right')
  table["logrho"] = table["logrho"][:irmax]
  for k in table.keys():
    if len(table[k].shape) == 3:
      table[k] = table[k][:,:,:irmax]

# Convert/rename variables
table["density"]        = 10.0**table["logrho"]
table["depsdT_rho"]     =       table["dedt"]
table["dpdeps_rhoye"]   =       table["dpderho"]
table["dpdrho_epsye"]   =       table["dpdrhoe"]
table["internalEnergy"] = 10.0**table["logenergy"]
table["pressure"]       = 10.0**table["logpress"]
table["temperature"]    = 10.0**table["logtemp"]
table["version"]        = np.array([1])
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Get rid of the energy shift
# -----------------------------------------------------------------------------
if args.mass_factor is None:
    table["mass_factor"] = ut.UAMU_MEV*(1.0 - table["energy_shift"]/ut.C_CGS**2)
else:
    table["mass_factor"] = args.mass_factor
table["mass_factor"] = np.array(table["mass_factor"])

table["internalEnergy"] -= table["energy_shift"]
table["energy_shift"][0] = 0.0

# energy density in units of c^2 g
rho3d = table["density"].reshape((1, 1, table["density"].shape[0]))
energy = rho3d * (1 + table["internalEnergy"]/(ut.C_CGS**2))

# number density in fm^{-3}
ndens = table["density"]/(ut.UAMU_CGS/(ut.FM_CGS**3))
mass_factor_cgs = (table["mass_factor"]*ut.MEV_CGS)/(ut.C_CGS**2)

# new density definition
table["density"] = ndens*(mass_factor_cgs/ut.FM_CGS**3)
table["logrho"] = np.log10(table["density"])

table["internalEnergy"] = ut.C_CGS**2 * (energy/table["density"] - 1.0)
assert np.all(table["internalEnergy"] > 0)
table["logenergy"] = np.log10(table["internalEnergy"])
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Compute relativistic sound speed
# -----------------------------------------------------------------------------
enthalpy = ut.C_CGS**2 + table["internalEnergy"] + table["pressure"]/rho3d
table["cs2"] /= enthalpy
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Cleanup table
# -----------------------------------------------------------------------------
# Clean sound speed
table["cs2"][np.isnan(table["cs2"])] = 1.0
table["cs2"][table["cs2"] < 0]       = 0.0
table["cs2"][table["cs2"] > 1]       = 1.0
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Output table
# -----------------------------------------------------------------------------
tstamp = datetime.date.today().strftime("%d-%b-%Y")

ofname = args.name + "_hydro_" + tstamp + ".h5"
dfile = h5py.File(ofname, "w")
for var in tabentries.hydro:
    if table[var].dtype == np.float64:
        assert np.all(np.isfinite(table[var]))
    dfile.create_dataset(var, data=table[var])
dfile.close()

ofname = args.name + "_weak_" + tstamp + ".h5"
dfile = h5py.File(ofname, "w")
for var in tabentries.weak:
    if table[var].dtype == np.float64:
        assert np.all(np.isfinite(table[var]))
    dfile.create_dataset(var, data=table[var])
dfile.close()

ofname = args.name + "_comp_" + tstamp + ".h5"
dfile = h5py.File(ofname, "w")
for var in tabentries.comp:
    if var in table.keys():
        if table[var].dtype == np.float64:
            assert np.all(np.isfinite(table[var]))
        dfile.create_dataset(var, data=table[var])
dfile.close()
# -----------------------------------------------------------------------------
