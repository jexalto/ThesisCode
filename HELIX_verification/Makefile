nocorrection:
    @if [ ! -f "helix/src/aerodynamics/correction-models/mod_prandtl.f90" ]; then \
    echo "------------- NO correction found -------------";\
    (cd HELIX_verification && python3 HELIX_Script_nocorr.py);\
    else \
    (cd helix && git checkout master);\
    (cd helix && make clean);\
    (cd helix && make HELIX_PYTHON);\
    (cd HELIX_verification && python3 HELIX_Script_nocorr.py);\
    fi;

correction:
    @if [ ! -f "helix/src/aerodynamics/correction-models/mod_prandtl.f90" ]; then \
    echo "------------- NO correction found -------------";\
    (cd helix && git checkout CorrectionFactorPrandtl);\
    (cd helix && make clean);\
    (cd helix && make HELIX_PYTHON);\
    (cd HELIX_verification && python3 HELIX_Script_corr.py);\
    else \
    (cd HELIX_verification && python3 HELIX_Script_corr.py);\
    fi;

total:
    @if [ ! -f "helix/src/aerodynamics/correction-models/mod_prandtl.f90" ]; then \
    (cd HELIX_verification && python3 HELIX_Script_nocorr.py);\
    echo " ------------- NO correction complete -------------";\
    (cd helix && git checkout CorrectionFactorPrandtl);\
    (cd helix && make clean);\
    (cd helix && make HELIX_PYTHON);\
    (cd HELIX_verification && python3 HELIX_Script_corr.py);\
    else \
    (cd HELIX_verification && python3 HELIX_Script_corr.py);\
    echo " ------------- Correction complete -------------";\
    (cd helix && git checkout master );\
    (cd helix && make clean);\
    (cd helix && make HELIX_PYTHON);\
    (cd HELIX_verification && python3 HELIX_Script_nocorr.py);\
    fi;

troubleshoot:
    (cd helix && make HELIX_PYTHON);\
    (cd HELIX_verification && python3 HELIX_Script_corr.py);\