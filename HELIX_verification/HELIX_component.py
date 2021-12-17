# Helix Libraries
import helix.model

import helix.parameters.simparam_def as simparam_def

import helix.geometry.geometry_def as geo_def
import helix.geometry.geometry_def_parametric as geo_def_parametric

import helix.references.references_def as ref_def

import helix_pyf90.mod_solver_initialize as f90_solver_initialize
import helix_pyf90.mod_aerodynamic_coefficients as f90_aerodynamic_coefficients
# External
import numpy as np

###############################################################################
## Main - Run Script
###############################################################################
def run_helix(rpm, v_inf):
    # ---------------------- Assign Simulation Parameters ------------------- #
    simparam = set_simparam(v_inf)

    # -------------------------- Set Reference Frames ----------------------- #
    references_def = set_references()

    # ----------------------------- Build Vehicle --------------------------- #
    geometry_def = geometry_definition(rpm)

    # ------------------------------- Run HELIX ----------------------------- #
    model = helix.model.Model()
    model.simparam_def = simparam
    model.geometry_def = geometry_def
    model.references_def = references_def
    model.initialize_run()
    model.run()
    return model.f_geometry.rotor[0].thrust_mag[0]


###############################################################################
## Set Simulation Parameters
###############################################################################
def set_simparam(v_inf):
    simparam = simparam_def.t_simparam_def()
    simparam.basename = "MR8x45 Rotor"

    simparam.dt = 0.5
    simparam.t_start = 0.0
    simparam.t_end = 1.0

    simparam.nt_rev = 30

    # Flipping freestream to make Z-axis of rotation
    simparam.v_inf = v_inf
    simparam.rho_inf = 1.25

    return simparam


###############################################################################
## Generate Vehicle
###############################################################################
def geometry_definition(rpm):
    # Initialize Geometry Component Definitions Holder
    geometry_def = geo_def.t_geometry_def()

    # ==========================================================================
    # Define Rotor 1
    # ==========================================================================
    # ---------------------------- Blade Parameters -------------------------- #
    rotor1 = geo_def_parametric.t_geometry_def_parametric()
    rotor1.CompName = "Rotor1"
    rotor1.CompType = "rotor"
    rotor1.RefName = "Hub"

    # Reference Parameters
    N_span = 35
    rotor1.ref_point = np.array([0.0, 0.02023364, 0.0])
    rotor1.ref_chord_frac = 0.5

    # Symmetry Parameters
    rotor1.symmetry = False
    rotor1.symmetry_point = np.array([0.0, 0.0, 0.0])
    rotor1.symmetry_normal = np.array([0.0, 1.0, 0.0])

    # Mirror Parameters
    rotor1.mirror = False
    rotor1.mirror_point = np.array([0.0, 0.0, 0.0])
    rotor1.mirror_normal = np.array([0.0, 1.0, 0.0])

    # Initialize Rotor and Allocate Arrays
    rotor1.initialize_parametric_geometry_definition(N_span)

    rotor1.multiple = True
    rotor1.multiplicity = {
        "mult_type": "rotor",
        "n_blades": 2,
        "rot_axis": np.array([0.0, 0.0, 1.0]),
        "rot_rate": rpm / 60.0 * 2.0 * np.pi,
        "psi_0": 0.0,
        "hub_offset": 0.0,
        "n_dofs": 0,
    }

    # ------------------------ Blade Section Definition ---------------------- #
    # Chord 1 ------------------
    rotor1.sec[0].chord = 0.01761998
    rotor1.sec[0].twist = 40.6848
    rotor1.sec[0].thick = 0.00216916
    rotor1.sec[0].alpha_0 = 0.3595378
    rotor1.sec[0].alpha_L0 = -0.03490658503
    rotor1.sec[0].Cl_alpha = 5.3407
    rotor1.sec[0].M = 50.0

    # Span 1  ------------------
    rotor1.span[0].span = 0.00126238
    rotor1.span[0].sweep = 0.0
    rotor1.span[0].dihed = 0.0
    rotor1.span[0].N_elem_span = 1
    rotor1.span[0].span_type = 1

    # Chord 2 ------------------
    rotor1.sec[1].chord = 0.0183515
    rotor1.sec[1].twist = 39.8284
    rotor1.sec[1].thick = 0.00211836
    rotor1.sec[1].alpha_0 = 0.3595378
    rotor1.sec[1].alpha_L0 = -0.03490658503
    rotor1.sec[1].Cl_alpha = 5.3407
    rotor1.sec[1].M = 50.0

    # Span 2  ------------------
    rotor1.span[1].span = 0.00126238
    rotor1.span[1].sweep = 0.0
    rotor1.span[1].dihed = 0.0
    rotor1.span[1].N_elem_span = 1
    rotor1.span[1].span_type = 1

    # Chord 3 ------------------
    rotor1.sec[2].chord = 0.01902714
    rotor1.sec[2].twist = 38.6086
    rotor1.sec[2].thick = 0.00207264
    rotor1.sec[2].alpha_0 = 0.3595378
    rotor1.sec[2].alpha_L0 = -0.03490658503
    rotor1.sec[2].Cl_alpha = 5.3407
    rotor1.sec[2].M = 50.0

    # Span 3  ------------------
    rotor1.span[2].span = 0.00126238
    rotor1.span[2].sweep = 0.0
    rotor1.span[2].dihed = 0.0
    rotor1.span[2].N_elem_span = 1
    rotor1.span[2].span_type = 1

    # Chord 4 ------------------
    rotor1.sec[3].chord = 0.01964182
    rotor1.sec[3].twist = 37.1367
    rotor1.sec[3].thick = 0.00202946
    rotor1.sec[3].alpha_0 = 0.3595378
    rotor1.sec[3].alpha_L0 = -0.03490658503
    rotor1.sec[3].Cl_alpha = 5.3407
    rotor1.sec[3].M = 50.0

    # Span 4  ------------------
    rotor1.span[3].span = 0.00126238
    rotor1.span[3].sweep = 0.0
    rotor1.span[3].dihed = 0.0
    rotor1.span[3].N_elem_span = 1
    rotor1.span[3].span_type = 1

    # Chord 5 ------------------
    rotor1.sec[4].chord = 0.02020062
    rotor1.sec[4].twist = 35.7341
    rotor1.sec[4].thick = 0.00199136
    rotor1.sec[4].alpha_0 = 0.3595378
    rotor1.sec[4].alpha_L0 = -0.03490658503
    rotor1.sec[4].Cl_alpha = 5.3407
    rotor1.sec[4].M = 50.0

    # Span 5  ------------------
    rotor1.span[4].span = 0.00126492
    rotor1.span[4].sweep = 0.0
    rotor1.span[4].dihed = 0.0
    rotor1.span[4].N_elem_span = 1
    rotor1.span[4].span_type = 1

    # Chord 6 ------------------
    rotor1.sec[5].chord = 0.02070354
    rotor1.sec[5].twist = 34.4208
    rotor1.sec[5].thick = 0.00196088
    rotor1.sec[5].alpha_0 = 0.3595378
    rotor1.sec[5].alpha_L0 = -0.03490658503
    rotor1.sec[5].Cl_alpha = 5.3407
    rotor1.sec[5].M = 50.0

    # Span 6  ------------------
    rotor1.span[5].span = 0.001778
    rotor1.span[5].sweep = 0.0
    rotor1.span[5].dihed = 0.0
    rotor1.span[5].N_elem_span = 1
    rotor1.span[5].span_type = 1

    # Chord 7 ------------------
    rotor1.sec[6].chord = 0.02131314
    rotor1.sec[6].twist = 32.7087
    rotor1.sec[6].thick = 0.0019304
    rotor1.sec[6].alpha_0 = 0.3595378
    rotor1.sec[6].alpha_L0 = -0.03490658503
    rotor1.sec[6].Cl_alpha = 5.3407
    rotor1.sec[6].M = 50.0

    # Span 7  ------------------
    rotor1.span[6].span = 0.00248666
    rotor1.span[6].sweep = 0.0
    rotor1.span[6].dihed = 0.0
    rotor1.span[6].N_elem_span = 1
    rotor1.span[6].span_type = 1

    # Chord 8 ------------------
    rotor1.sec[7].chord = 0.02197608
    rotor1.sec[7].twist = 30.5564
    rotor1.sec[7].thick = 0.00189992
    rotor1.sec[7].alpha_0 = 0.3595378
    rotor1.sec[7].alpha_L0 = -0.03490658503
    rotor1.sec[7].Cl_alpha = 5.3407
    rotor1.sec[7].M = 50.0

    # Span 8  ------------------
    rotor1.span[7].span = 0.00250698
    rotor1.span[7].sweep = 0.0
    rotor1.span[7].dihed = 0.0
    rotor1.span[7].N_elem_span = 1
    rotor1.span[7].span_type = 1

    # Chord 9 ------------------
    rotor1.sec[8].chord = 0.02242566
    rotor1.sec[8].twist = 28.6333
    rotor1.sec[8].thick = 0.00188214
    rotor1.sec[8].alpha_0 = 0.3595378
    rotor1.sec[8].alpha_L0 = -0.03490658503
    rotor1.sec[8].Cl_alpha = 5.3407
    rotor1.sec[8].M = 50.0

    # Span 9  ------------------
    rotor1.span[8].span = 0.00250444
    rotor1.span[8].sweep = 0.0
    rotor1.span[8].dihed = 0.0
    rotor1.span[8].N_elem_span = 1
    rotor1.span[8].span_type = 1

    # Chord 10 -----------------
    rotor1.sec[9].chord = 0.02265172
    rotor1.sec[9].twist = 26.9206
    rotor1.sec[9].thick = 0.00186182
    rotor1.sec[9].alpha_0 = 0.3595378
    rotor1.sec[9].alpha_L0 = -0.03490658503
    rotor1.sec[9].Cl_alpha = 5.3407
    rotor1.sec[9].M = 50.0

    # Span 10  ------------------
    rotor1.span[9].span = 0.00250698
    rotor1.span[9].sweep = 0.0
    rotor1.span[9].dihed = 0.0
    rotor1.span[9].N_elem_span = 1
    rotor1.span[9].span_type = 1

    # Chord 11 -----------------
    rotor1.sec[10].chord = 0.02265934
    rotor1.sec[10].twist = 25.3882
    rotor1.sec[10].thick = 0.00182626
    rotor1.sec[10].alpha_0 = 0.3595378
    rotor1.sec[10].alpha_L0 = -0.03490658503
    rotor1.sec[10].Cl_alpha = 5.3407
    rotor1.sec[10].M = 50.0

    # Span 11  ------------------
    rotor1.span[10].span = 0.00250698
    rotor1.span[10].sweep = 0.0
    rotor1.span[10].dihed = 0.0
    rotor1.span[10].N_elem_span = 1
    rotor1.span[10].span_type = 1

    # Chord 12 -----------------
    rotor1.sec[11].chord = 0.02248662
    rotor1.sec[11].twist = 24.0111
    rotor1.sec[11].thick = 0.00178054
    rotor1.sec[11].alpha_0 = 0.3595378
    rotor1.sec[11].alpha_L0 = -0.03490658503
    rotor1.sec[11].Cl_alpha = 5.3407
    rotor1.sec[11].M = 50.0

    # Span 12  ------------------
    rotor1.span[11].span = 0.00250444
    rotor1.span[11].sweep = 0.0
    rotor1.span[11].dihed = 0.0
    rotor1.span[11].N_elem_span = 1
    rotor1.span[11].span_type = 1

    # Chord 13 -----------------
    rotor1.sec[12].chord = 0.02225548
    rotor1.sec[12].twist = 22.7681
    rotor1.sec[12].thick = 0.00173482
    rotor1.sec[12].alpha_0 = 0.3595378
    rotor1.sec[12].alpha_L0 = -0.03490658503
    rotor1.sec[12].Cl_alpha = 5.3407
    rotor1.sec[12].M = 50.0

    # Span 13  ------------------
    rotor1.span[12].span = 0.00250698
    rotor1.span[12].sweep = 0.0
    rotor1.span[12].dihed = 0.0
    rotor1.span[12].N_elem_span = 1
    rotor1.span[12].span_type = 1

    # Chord 14 -----------------
    rotor1.sec[13].chord = 0.021971
    rotor1.sec[13].twist = 21.6415
    rotor1.sec[13].thick = 0.00168656
    rotor1.sec[13].alpha_0 = 0.3595378
    rotor1.sec[13].alpha_L0 = -0.03490658503
    rotor1.sec[13].Cl_alpha = 5.3407
    rotor1.sec[13].M = 50.0

    # Span 14  ------------------
    rotor1.span[13].span = 0.00250444
    rotor1.span[13].sweep = 0.0
    rotor1.span[13].dihed = 0.0
    rotor1.span[13].N_elem_span = 1
    rotor1.span[13].span_type = 1

    # Chord 15 -----------------
    rotor1.sec[14].chord = 0.0216408
    rotor1.sec[14].twist = 20.6165
    rotor1.sec[14].thick = 0.00164084
    rotor1.sec[14].alpha_0 = 0.3595378
    rotor1.sec[14].alpha_L0 = -0.03490658503
    rotor1.sec[14].Cl_alpha = 5.3407
    rotor1.sec[14].M = 50.0

    # Span 15  ------------------
    rotor1.span[14].span = 0.00250698
    rotor1.span[14].sweep = 0.0
    rotor1.span[14].dihed = 0.0
    rotor1.span[14].N_elem_span = 1
    rotor1.span[14].span_type = 1

    # Chord 16 -----------------
    rotor1.sec[15].chord = 0.02126488
    rotor1.sec[15].twist = 19.6805
    rotor1.sec[15].thick = 0.00159258
    rotor1.sec[15].alpha_0 = 0.3595378
    rotor1.sec[15].alpha_L0 = -0.03490658503
    rotor1.sec[15].Cl_alpha = 5.3407
    rotor1.sec[15].M = 50.0

    # Span 16 ------------------
    rotor1.span[15].span = 0.00250698
    rotor1.span[15].sweep = 0.0
    rotor1.span[15].dihed = 0.0
    rotor1.span[15].N_elem_span = 1
    rotor1.span[15].span_type = 1

    # Chord 17 -----------------
    rotor1.sec[16].chord = 0.02084324
    rotor1.sec[16].twist = 18.8228
    rotor1.sec[16].thick = 0.00154686
    rotor1.sec[16].alpha_0 = 0.3595378
    rotor1.sec[16].alpha_L0 = -0.03490658503
    rotor1.sec[16].Cl_alpha = 5.3407
    rotor1.sec[16].M = 50.0

    # Span 17 ------------------
    rotor1.span[16].span = 0.00250444
    rotor1.span[16].sweep = 0.0
    rotor1.span[16].dihed = 0.0
    rotor1.span[16].N_elem_span = 1
    rotor1.span[16].span_type = 1

    # Chord 18 -----------------
    rotor1.sec[17].chord = 0.02037842
    rotor1.sec[17].twist = 18.0344
    rotor1.sec[17].thick = 0.0014986
    rotor1.sec[17].alpha_0 = 0.3595378
    rotor1.sec[17].alpha_L0 = -0.03490658503
    rotor1.sec[17].Cl_alpha = 5.3407
    rotor1.sec[17].M = 50.0

    # Span 18  ------------------
    rotor1.span[17].span = 0.00250698
    rotor1.span[17].sweep = 0.0
    rotor1.span[17].dihed = 0.0
    rotor1.span[17].N_elem_span = 1
    rotor1.span[17].span_type = 1

    # Chord 19 -----------------
    rotor1.sec[18].chord = 0.01987296
    rotor1.sec[18].twist = 17.3074
    rotor1.sec[18].thick = 0.00145034
    rotor1.sec[18].alpha_0 = 0.3595378
    rotor1.sec[18].alpha_L0 = -0.03490658503
    rotor1.sec[18].Cl_alpha = 5.3407
    rotor1.sec[18].M = 50.0

    # Span 19  ------------------
    rotor1.span[18].span = 0.00250698
    rotor1.span[18].sweep = 0.0
    rotor1.span[18].dihed = 0.0
    rotor1.span[18].N_elem_span = 1
    rotor1.span[18].span_type = 1

    # Chord 20 -----------------
    rotor1.sec[19].chord = 0.01932686
    rotor1.sec[19].twist = 16.6352
    rotor1.sec[19].thick = 0.00140208
    rotor1.sec[19].alpha_0 = 0.3595378
    rotor1.sec[19].alpha_L0 = -0.03490658503
    rotor1.sec[19].Cl_alpha = 5.3407
    rotor1.sec[19].M = 50.0

    # Span 20  ------------------
    rotor1.span[19].span = 0.00250444
    rotor1.span[19].sweep = 0.0
    rotor1.span[19].dihed = 0.0
    rotor1.span[19].N_elem_span = 1
    rotor1.span[19].span_type = 1

    # Chord 21 -----------------
    rotor1.sec[20].chord = 0.01874266
    rotor1.sec[20].twist = 16.0119
    rotor1.sec[20].thick = 0.00135382
    rotor1.sec[20].alpha_0 = 0.3595378
    rotor1.sec[20].alpha_L0 = -0.03490658503
    rotor1.sec[20].Cl_alpha = 5.3407
    rotor1.sec[20].M = 50.0

    # Span 21  ------------------
    rotor1.span[20].span = 0.00250698
    rotor1.span[20].sweep = 0.0
    rotor1.span[20].dihed = 0.0
    rotor1.span[20].N_elem_span = 1
    rotor1.span[20].span_type = 1

    # Chord 22 -----------------
    rotor1.sec[21].chord = 0.0181229
    rotor1.sec[21].twist = 15.4326
    rotor1.sec[21].thick = 0.00130302
    rotor1.sec[21].alpha_0 = 0.3595378
    rotor1.sec[21].alpha_L0 = -0.03490658503
    rotor1.sec[21].Cl_alpha = 5.3407
    rotor1.sec[21].M = 50.0

    # Span 22  ------------------
    rotor1.span[21].span = 0.00250444
    rotor1.span[21].sweep = 0.0
    rotor1.span[21].dihed = 0.0
    rotor1.span[21].N_elem_span = 1
    rotor1.span[21].span_type = 1

    # Chord 23 -----------------
    rotor1.sec[22].chord = 0.01746504
    rotor1.sec[22].twist = 14.8928
    rotor1.sec[22].thick = 0.00125222
    rotor1.sec[22].alpha_0 = 0.3595378
    rotor1.sec[22].alpha_L0 = -0.03490658503
    rotor1.sec[22].Cl_alpha = 5.3407
    rotor1.sec[22].M = 50.0

    # Span 23  ------------------
    rotor1.span[22].span = 0.00250698
    rotor1.span[22].sweep = 0.0
    rotor1.span[22].dihed = 0.0
    rotor1.span[22].N_elem_span = 1
    rotor1.span[22].span_type = 1

    # Chord 24 -----------------
    rotor1.sec[23].chord = 0.0167767
    rotor1.sec[23].twist = 14.3887
    rotor1.sec[23].thick = 0.00120142
    rotor1.sec[23].alpha_0 = 0.3595378
    rotor1.sec[23].alpha_L0 = -0.03490658503
    rotor1.sec[23].Cl_alpha = 5.3407
    rotor1.sec[23].M = 50.0

    # Span 24  ------------------
    rotor1.span[23].span = 0.00250698
    rotor1.span[23].sweep = 0.0
    rotor1.span[23].dihed = 0.0
    rotor1.span[23].N_elem_span = 1
    rotor1.span[23].span_type = 1

    # Chord 25 -----------------
    rotor1.sec[24].chord = 0.01605788
    rotor1.sec[24].twist = 13.9169
    rotor1.sec[24].thick = 0.00115062
    rotor1.sec[24].alpha_0 = 0.3595378
    rotor1.sec[24].alpha_L0 = -0.03490658503
    rotor1.sec[24].Cl_alpha = 5.3407
    rotor1.sec[24].M = 50.0

    # Span 25  ------------------
    rotor1.span[24].span = 0.00250444
    rotor1.span[24].sweep = 0.0
    rotor1.span[24].dihed = 0.0
    rotor1.span[24].N_elem_span = 1
    rotor1.span[24].span_type = 1

    # Chord 26 -----------------
    rotor1.sec[25].chord = 0.01530604
    rotor1.sec[25].twist = 13.4746
    rotor1.sec[25].thick = 0.00109728
    rotor1.sec[25].alpha_0 = 0.3595378
    rotor1.sec[25].alpha_L0 = -0.03490658503
    rotor1.sec[25].Cl_alpha = 5.3407
    rotor1.sec[25].M = 50.0

    # Span 26  ------------------
    rotor1.span[25].span = 0.00250698
    rotor1.span[25].sweep = 0.0
    rotor1.span[25].dihed = 0.0
    rotor1.span[25].N_elem_span = 1
    rotor1.span[25].span_type = 1

    # Chord 27 -----------------
    rotor1.sec[26].chord = 0.0145288
    rotor1.sec[26].twist = 13.059
    rotor1.sec[26].thick = 0.0010414
    rotor1.sec[26].alpha_0 = 0.3595378
    rotor1.sec[26].alpha_L0 = -0.03490658503
    rotor1.sec[26].Cl_alpha = 5.3407
    rotor1.sec[26].M = 50.0

    # Span 27  ------------------
    rotor1.span[26].span = 0.00250444
    rotor1.span[26].sweep = 0.0
    rotor1.span[26].dihed = 0.0
    rotor1.span[26].N_elem_span = 1
    rotor1.span[26].span_type = 1

    # Chord 28 -----------------
    rotor1.sec[27].chord = 0.01372362
    rotor1.sec[27].twist = 12.6679
    rotor1.sec[27].thick = 0.00098552
    rotor1.sec[27].alpha_0 = 0.3595378
    rotor1.sec[27].alpha_L0 = -0.03490658503
    rotor1.sec[27].Cl_alpha = 5.3407
    rotor1.sec[27].M = 50.0

    # Span 28  ------------------
    rotor1.span[27].span = 0.00250698
    rotor1.span[27].sweep = 0.0
    rotor1.span[27].dihed = 0.0
    rotor1.span[27].N_elem_span = 1
    rotor1.span[27].span_type = 1

    # Chord 29 -----------------
    rotor1.sec[28].chord = 0.01289304
    rotor1.sec[28].twist = 12.2992
    rotor1.sec[28].thick = 0.0009271
    rotor1.sec[28].alpha_0 = 0.3595378
    rotor1.sec[28].alpha_L0 = -0.03490658503
    rotor1.sec[28].Cl_alpha = 5.3407
    rotor1.sec[28].M = 50.0

    # Span 29  ------------------
    rotor1.span[28].span = 0.00250698
    rotor1.span[28].sweep = 0.0
    rotor1.span[28].dihed = 0.0
    rotor1.span[28].N_elem_span = 1
    rotor1.span[28].span_type = 1

    # Chord 30 -----------------
    rotor1.sec[29].chord = 0.0120396
    rotor1.sec[29].twist = 11.951
    rotor1.sec[29].thick = 0.00086868
    rotor1.sec[29].alpha_0 = 0.3595378
    rotor1.sec[29].alpha_L0 = -0.03490658503
    rotor1.sec[29].Cl_alpha = 5.3407
    rotor1.sec[29].M = 50.0

    # Span 30  ------------------
    rotor1.span[29].span = 0.00250444
    rotor1.span[29].sweep = 0.0
    rotor1.span[29].dihed = 0.0
    rotor1.span[29].N_elem_span = 1
    rotor1.span[29].span_type = 1

    # Chord 31 -----------------
    rotor1.sec[30].chord = 0.0111633
    rotor1.sec[30].twist = 11.6218
    rotor1.sec[30].thick = 0.00080772
    rotor1.sec[30].alpha_0 = 0.3595378
    rotor1.sec[30].alpha_L0 = -0.03490658503
    rotor1.sec[30].Cl_alpha = 5.3407
    rotor1.sec[30].M = 50.0

    # Span 31  ------------------
    rotor1.span[30].span = 0.00250698
    rotor1.span[30].sweep = 0.0
    rotor1.span[30].dihed = 0.0
    rotor1.span[30].N_elem_span = 1
    rotor1.span[30].span_type = 1

    # Chord 32 -----------------
    rotor1.sec[31].chord = 0.01026668
    rotor1.sec[31].twist = 11.3099
    rotor1.sec[31].thick = 0.00074676
    rotor1.sec[31].alpha_0 = 0.3595378
    rotor1.sec[31].alpha_L0 = -0.03490658503
    rotor1.sec[31].Cl_alpha = 5.3407
    rotor1.sec[31].M = 50.0

    # Span 32  ------------------
    rotor1.span[31].span = 0.00250444
    rotor1.span[31].sweep = 0.0
    rotor1.span[31].dihed = 0.0
    rotor1.span[31].N_elem_span = 1
    rotor1.span[31].span_type = 1

    # Chord 33 -----------------
    rotor1.sec[32].chord = 0.00935228
    rotor1.sec[32].twist = 11.0142
    rotor1.sec[32].thick = 0.00068072
    rotor1.sec[32].alpha_0 = 0.3595378
    rotor1.sec[32].alpha_L0 = -0.03490658503
    rotor1.sec[32].Cl_alpha = 5.3407
    rotor1.sec[32].M = 50.0

    # Span 33  ------------------
    rotor1.span[32].span = 0.00250698
    rotor1.span[32].sweep = 0.0
    rotor1.span[32].dihed = 0.0
    rotor1.span[32].N_elem_span = 1
    rotor1.span[32].span_type = 1

    # Chord 34 -----------------
    rotor1.sec[33].chord = 0.00842264
    rotor1.sec[33].twist = 10.7334
    rotor1.sec[33].thick = 0.00061722
    rotor1.sec[33].alpha_0 = 0.3595378
    rotor1.sec[33].alpha_L0 = -0.03490658503
    rotor1.sec[33].Cl_alpha = 5.3407
    rotor1.sec[33].M = 50.0

    # Span 34  ------------------
    rotor1.span[33].span = 0.0023368
    rotor1.span[33].sweep = 0.0
    rotor1.span[33].dihed = 0.0
    rotor1.span[33].N_elem_span = 1
    rotor1.span[33].span_type = 1

    # Chord 35 -----------------
    rotor1.sec[34].chord = 0.00705866
    rotor1.sec[34].twist = 10.4841
    rotor1.sec[34].thick = 0.00051816
    rotor1.sec[34].alpha_0 = 0.3595378
    rotor1.sec[34].alpha_L0 = -0.03490658503
    rotor1.sec[34].Cl_alpha = 5.3407
    rotor1.sec[34].M = 50.0

    # Span 35  ------------------
    rotor1.span[34].span = 0.002286
    rotor1.span[34].sweep = 0.0
    rotor1.span[34].dihed = 0.0
    rotor1.span[34].N_elem_span = 1
    rotor1.span[34].span_type = 1

    # Chord 36 -----------------
    rotor1.sec[35].chord = 0.00395478
    rotor1.sec[35].twist = 10.2509
    rotor1.sec[35].thick = 0.0002921
    rotor1.sec[35].alpha_0 = 0.3595378
    rotor1.sec[35].alpha_L0 = -0.03490658503
    rotor1.sec[35].Cl_alpha = 5.3407
    rotor1.sec[35].M = 50.0

    # Append To Vehicle
    geometry_def.append_component(rotor1)

    return geometry_def


# ###############################################################################
# ## Generate Vehicle
# ###############################################################################
def set_references():
    # Initialize Reference Frame Defintions Holder
    references_def = ref_def.t_references_def()

    # ==========================================================================
    # Hub Frame
    # ==========================================================================
    Hub = ref_def.t_frame_def()
    Hub.Name = "Hub"
    Hub.Parent = "Root"
    Hub.origin = np.array([0.0, 0.0, 0.0])
    Hub.orientation = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    Hub.moving = False

    # Append to References
    references_def.append_frame(Hub)

    return references_def