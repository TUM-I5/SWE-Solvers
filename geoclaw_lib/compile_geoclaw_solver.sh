#!/bin/bash

GEOCLAWSRCPATH=../src/solver/geoclaw

ifort -w -g -c $GEOCLAWSRCPATH/riemannsolvers.f $GEOCLAWSRCPATH/c_bind_riemannsolvers.f90