clear
% Load 5-bus system without line 1-4
bus5 = load('bus5.mat');

%% ROM for bus 1

inputs = [1 2 3 4];
outputs = inputs;
sp = SingularPerturbation(bus5.sys(inputs, outputs));

% 1:20  = GFLI #1 + GFLI #2
% 41:46 = inf bus + shunts 1 + shunt 2
% 51:54 = shunt 5 + line12
% 62:66 = line15
bus1_states = [1:20 41:46 51:54 62:66];
rstates = setdiff(1:66, bus1_states);
"Bus 1 ROM with "+length(bus1_states)+" states"
rom1 = getrom(sp, rstates, "first");
error(rom1, Inf)
sigmaplot(rom1)
%% ROM for bus 3

inputs = [5 6];
outputs = inputs;
sp = SingularPerturbation(bus5.sys(inputs, outputs));

% 21:30  = GFLI #3
% 45:50 = shunt 2/3/4
% 55:58 = line23 + line 34
bus3_states = [21:30 45:50 55:58];
rstates = setdiff(1:66, bus3_states);
"Bus 3 ROM with "+length(bus3_states)+" states"
rom3 = getrom(sp, rstates, "zero");
error(rom3, Inf)
bodeplot(rom3)
sigmaplot(rom3)

%% ROM for bus 4

inputs = [7 8 9];
outputs = [7 8];
sp = SingularPerturbation(bus5.sys(inputs, outputs));

% 31:42  = GFLI #3 + inf bus
% 47:52 = shunt 3/4/5
% 57:62 = line34 + line 45
bus4_states = [31:42 47:52 55:62];
rstates = setdiff(1:66, bus4_states);
"Bus 4 ROM with "+length(bus4_states)+" states"
rom4 = getrom(sp, rstates, "zero");
error(rom4, Inf)
bodeplot(rom4)
%sigmaplot(rom4)

%% TODO: 
% - Check eigenvalues of ROM
% - Iteratively reducing should not affect ROM, check this
% - Perform reduction for the 5-bus system
%   * Remove line from the 5-bus for more radial structure
% - Hermite interpolate
% What frequencies are relevant in dq-frame