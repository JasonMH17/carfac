% // clang-format off
% Copyright 2012 The CARFAC Authors. All Rights Reserved.
% Author: Richard F. Lyon
%
% This file is part of an implementation of Lyon's cochlear model:
% "Cascade of Asymmetric Resonators with Fast-Acting Compression"
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function CF = CARFAC_Init(CF)
% function CF = CARFAC_Init(CF)
%
% Initialize state for one or more ears of CF.
% This allocates and zeros all the state vector storage in the CF struct.
% The CAR state starts at rest for the current OHC_health, so after changing
% OHC_health call this again and there is no start-up transient.

n_ears = CF.n_ears;

for ear = 1:n_ears
  % for now there's only one coeffs, not one per ear
  CF.ears(ear).CAR_state = CAR_Init_State(CF.ears(ear).CAR_coeffs);
  CF.ears(ear).IHC_state = IHC_Init_State(CF.ears(ear).IHC_coeffs);
  CF.ears(ear).AGC_state = AGC_Init_State(CF.ears(ear).AGC_coeffs);
  if CF.do_syn
     CF.ears(ear).SYN_state = SYN_Init_State(CF.ears(ear).SYN_coeffs);
  end
end


function state = CAR_Init_State(coeffs)
% At rest the AGC memory is 0, so CARFAC_Close_AGC_Loop holds the undamping at
% OHC_health, with zB_memory at zr_coeffs .* OHC_health and g_memory at
% CARFAC_Stage_g of it.  Starting there, rather than at the healthy design's
% undamping of 1, avoids a start-up transient when OHC_health was changed
% after design.  CARFAC_Stage_g is the quadratic approximation that
% CARFAC_Close_AGC_Loop uses, so this is its rest for any OHC_health; with
% OHC_health of 1 it matches the design's g0_coeffs to rounding.
n_ch = coeffs.n_ch;
undamping = coeffs.OHC_health;
state = struct( ...
  'z1_memory', zeros(n_ch, 1), ...
  'z2_memory', zeros(n_ch, 1), ...
  'zA_memory', zeros(n_ch, 1), ...
  'zB_memory', coeffs.zr_coeffs .* undamping, ...
  'dzB_memory', zeros(n_ch, 1), ...
  'zY_memory', zeros(n_ch, 1), ...
  'g_memory', CARFAC_Stage_g(coeffs, undamping), ...
  'dg_memory', zeros(n_ch, 1), ...
  'ac_coupler', zeros(n_ch, 1) ...
  );


function state = AGC_Init_State(coeffs)
% 2025 new way, one struct instead of array of them.
state = struct( ...
  'AGC_memory', zeros(coeffs.n_ch, coeffs.n_AGC_stages), ...
   'decim_phase', zeros(1, coeffs.n_AGC_stages));  % small ints
if coeffs.simpler_decimating  % One decimation factor vs per stage.
  state.input_accum = zeros(coeffs.n_ch, 1);
else
  state.input_accum = zeros(coeffs.n_ch, coeffs.n_AGC_stages);
end


function state = IHC_Init_State(coeffs)
n_ch = coeffs.n_ch;
if coeffs.just_hwr
  state = struct('ihc_accum', zeros(n_ch, 1));
else
  if coeffs.one_cap
    state = struct( ...
      'ihc_accum', zeros(n_ch, 1), ...
      'cap_voltage', coeffs.rest_cap * ones(n_ch, 1), ...
      'lpf1_state', coeffs.rest_output * ones(n_ch, 1), ...
      'lpf2_state', coeffs.rest_output * ones(n_ch, 1) ...
      );
  else
    state = struct( ...
      'ihc_accum', zeros(n_ch, 1), ...
      'cap1_voltage', coeffs.rest_cap1 * ones(n_ch, 1), ...
      'cap2_voltage', coeffs.rest_cap2 * ones(n_ch, 1), ...
      'lpf1_state', coeffs.rest_output * ones(n_ch, 1) ...
      );
  end
end

function state = SYN_Init_State(coeffs)
n_ch = coeffs.n_ch;
n_cl = coeffs.n_classes;
state = struct( ...
  'reservoirs', ones(n_ch, 1) * coeffs.res_lpf_inits, ...  % 0 full, 1 empty.
  'lpf_state', ones(n_ch, 1) * coeffs.spont_p);
