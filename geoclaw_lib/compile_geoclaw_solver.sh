#!/bin/bash

GEOCLAWSRCPATH=../../geoclaw/src/2d/

ifort -w -g -c -fpp $GEOCLAWSRCPATH/riemannsolvers.f $GEOCLAWSRCPATH/c_bind_riemannsolvers.f90