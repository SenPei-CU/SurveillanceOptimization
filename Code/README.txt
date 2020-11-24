# SurveillanceOptimization
Code for Pei et al. Optimizing respiratory virus surveillance networks using uncertainty propagation

The code generates forecasting for a synthetic influenza outbreak in 35 US states.

# System requirements

Run the code in MATLAB (Mac or Windows). The code has been tested on MATLAB R2020a.

# Installation guide

Use standalone MATLAB. No need to install additional softwares.

# Demo

This code generates 1- to 4-week ahead ILI+ predictions starting from fweek.

For instance, to generate forecasts at week 10, call function "forecast(10)" in the MATLAB command window.

The code outputs a figure comparing forecasts and truths, as well as a mat file storing forecast results.

The expected run time for demo is within 1 minute on a typical desktop computer.

# Instruction for use

Users can replace the "Outbreakdata.mat" file with other outbreak data.

Surveillance networks can be modified by changing the variable "observed" in forecast.m.