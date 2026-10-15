#!/bin/bash
# Sourced by bin/pds-run-case.sbatch right before muscle_manager starts; every actor
# inherits it.
# my_contr writes its 14 voltages in coil order, the order the plasmaless actor reads
# (controllers/my_contr/muscle_NICE_input.m); its default, power-supply order, is NICE's.
export MY_CONTR_VOLTAGE_ORDER=coil
