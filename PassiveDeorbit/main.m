clear; clc; close all;

%% Disposal 
%% ATLANTE + 2xKIKI + 2xRANGOS disposal
Cd_Atlante = 2.2;          
mass_Atlante_kg = 4805.47;
mass_Kiki_kg = 606.58;
mass_Rangos_kg = 230;
area_Atlante_m2 = 100;  % this values are not the sum of the areas but are 
                           % constrained by the configuration, for now is
                           % set to solar panels' area

[disposalAtlante, burnAtlante, figAtlante] = ...
    PassiveDisposal(Cd_Atlante, mass_Atlante_kg+2*mass_Kiki_kg+2*mass_Rangos_kg, area_Atlante_m2);
%% TAXI disposal
Cd_Taxi = 2.2;            
mass_Taxi_kg = 1059.88;       
area_Taxi_m2 = 15;       %idk  

[disposalTaxi, burnTaxi, figTaxi] = ...
    PassiveDisposal(Cd_Taxi, mass_Taxi_kg, area_Taxi_m2);

%% Phasing with differential drag
% Function to be developed: PhasingDrag

% PURPOSE
% Assess whether differential drag can reduce the propulsive phasing cost
% within the available operational time.
%
% INPUTS
% - Chaser(Taxi): mass, Cd, minimum and maximum effective drag areas.
% - Target(Modules/Container): mass, Cd and effective drag area.
% - Initial orbital altitudes and initial along-track angular separation.
% - Maximum allowed phasing time and acceptable final angular separation.
% - Reference propulsive phasing DeltaV for the same final conditions.
%
% MODEL
% - Initially circular, coplanar orbits.
% - Propagate both semi-major axes under atmospheric drag.
% - Integrate the relative angular motion using the evolving mean motions.
% - Compare fixed-area and high/low-area chaser attitude strategies.
% - Repeat for the six reference altitudes and both atmospheric profiles.
% - Specify whether the target decays naturally or maintains its orbit.
%
% OUTPUTS
% - Table: phasing time, altitude loss, final orbital mismatch,
%   estimated correction DeltaV and net DeltaV saving.
% - Flag cases that exceed the allowed time or altitude-loss limits.
% - Plots: angular separation and orbital altitudes versus time.
%
% IMPORTANT
% Angular alignment alone does not complete rendezvous.
% Account for the final orbital mismatch before claiming a DeltaV saving.
% Area changes require an achievable attitude and compatible power conditions.

%% RAAN alignment through J2
% Function to be developed: RAANJ2
%
% PURPOSE
% Assess whether differential J2 nodal precession can reduce the propulsive
% cost of RAAN alignment within the available operational time.
%
% INPUTS
% - Chaser(Taxi) and target(Modules/Container): initial altitude, eccentricity, inclination and RAAN.
% - Candidate temporary drift-orbit altitudes for the chaser. (necessario se i drift
% delle RAAN dei componenti non producono (nelle orbite attuali) un
% deltaΩ=0 o che si avvicini a 0
% - Maximum allowed alignment time and final RAAN tolerance.
% - Reference propulsive alignment DeltaV for the same final conditions.
%
% MODEL
% - Compute the J2 nodal precession rates of both orbital planes.
% - Propagate their relative RAAN, accounting for the drift direction.
% - First assess alignment using the existing orbits without extra transfers.
% - Then assess temporary drift orbits, including entry and exit manoeuvres.
% - Repeat for the six reference orbits with their associated inclinations.
% - Initially assume constant orbital elements except RAAN.
% - Subsequently assess drag-induced altitude changes where significant.
%
% OUTPUTS
% - Table: relative precession rate, alignment time, drift-orbit altitude,
%   transfer DeltaV, residual correction DeltaV and net DeltaV saving.
% - Flag cases that exceed the allowed time or fail to reach alignment.
% - Plots: relative RAAN versus time and total DeltaV versus alignment time.
%
% IMPORTANT
% Both planes precess: alignment depends on their relative precession rate.
% Equal precession rates do not reduce the initial RAAN separation.
% RAAN alignment does not remove an inclination difference or solve phasing.
% Include drift-orbit transfers before claiming a DeltaV saving.
% Compare strategies at the same final orbit and target epoch.