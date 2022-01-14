This folder comprises the aero-propulsive model

runme.py is the main file of interest for the integration with the initiator. This is loaded to matlab as a module. This module is called in @AeroPropulsiveLoad and @AeroPropulsiveDeltas to calculate the aerodynamics of the system. It requires a matfile containing information about the design and flight conditions. This matfile is produced in the method 'getDesign' within the initiator modules and is stored as a temporary file. 

propellerperformance.py is only usd within the wing level study to estimate the required rpm to deliver a certain thrust.

BEM and VLM are the propeller and wing  models.
