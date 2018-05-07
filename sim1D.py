#!/usr/bin/python
# F. Massonnet, UCL, August 2017
#
# Set of functions intended to play with a simple sea-ice slab-ocean model.
# The (thermodynamic) sea ice model is a zero-layer (no thermal inertia) model
# without snow. The model is coupled to a mixed layer with uniform temperature.
# The atmospheric forcings are prescribed.
#
# This file consists of three parts:
# (1) A function to integrate the sea ice-ocean state by one time step, given a
#     forcing and a set of parameters
# (2) A function producing a default run. It calls (1) from a given initial 
#     state with default forcing and parameters.
# (3) A function to run (2) from outside Python (i.e. from shell)
#

# Import modules, etc.
# --------------------
import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime

# -----------------------------------------------------------------
# (1) The function to integrate the physical state by one time step
# -----------------------------------------------------------------

def sim1D(state_now, parameters, forcing, dt = 24.0 * 3600):
  """
  Function sim1D
     PURPOSE
       Takes the sea ice/ocean state now and returns the state at the next
       time step, given a set of parameters, a forcing, and a time step.
     INPUTS
       state_now : A Numpy array of dimension 4:
                   state_now[0] = Volume of ice per grid cell area [m3/m2]
                   state_now[1] = Area of ice per grid cell area   [m2/m2]
                   state_now[2] = Surface temperature of ice       [K]
                   state_now[3] = Ocean mixed layer temperature    [K]

       parameters : A Numpy array of dimension 11:
                    parameters[0]  = Stefan-Boltzman constant                  [W/m2/K4]
                    parameters[1]  = Freezing point of seawater                [K]
                    parameters[2]  = Specific seawater heat capacity           [K]
                    parameters[3]  = Seawater density                          [kg/m3]
                    parameters[4]  = Mixed layer depth                         [m]
                    parameters[5]  = Minimal sea ice thickness at creation     [m]
                    parameters[6]  = Volumetric latent heat of fusion of ice   [J/m3]
                    parameters[7]  = Sea ice heat conductivity                 [W/m/K]
                    parameters[8]  = Sea ice albedo                            [-]
                    parameters[9]  = Temperature of fusion of ice              [K]
                    parameters[10] = Albedo of ocean                           [-]

       forcing :    A Numpy array of dimension 2:
                    forcing[0] = Incoming shortwave radiation (below clouds)  [W/m2]
                    forcing[1] = Other fluxes (latent, sensible, long-wave)   [W/m2]
       
       dt :         Time step [s]

    OUTPUTS
      state_new : A Numpy array of dimension 4 by 1, giving the state updated
    REFERENCES
      The sea ice-ocean model was coded from scratch but largely follows 
        Fichefet, T. and Morales Maqueda, M. A. (1997), Sensitivity of a glocal sea ice model to the treatment 
        of the thermodynamics and dynamics, J. Geophys. Res., 102, C6, 12609-12646
  """

  # Get the parameters
  # Parameter                   Meaning                          Units             Recommended value   Source
  sigma     = parameters[0]   # Stefan-Boltzman constant.        W/m2/K4           5.67 * 10 ** (-8)   That's the way it is
  T_f       = parameters[1]   # Freezing point of seawater       K                 -1.8 + 273.15       
  cp_oce    = parameters[2]   # Specific seawater heat capacity  J/kg/K            4001.896            Fichefet and Morales Maqueda 1997           
  rho_oce   = parameters[3]   # Seawater density                 kg/m3             1024.458            F&MM 1997
  h_oce     = parameters[4]   # Mixed layer depth                m                 30.0                
  h_ice_min = parameters[5]   # Minimal sea ice thickness        m                 0.1                 F&MM 1997
  Lf_ice    = parameters[6]   # Latent heat of fusion of ice     J/m3              300.330 * 1e6       F&MM 1997
  k_ice     = parameters[7]   # Ice conductivity                 W/m/K             2.0344              F&MM 1997
  alb_ice   = parameters[8]   # Ice albedo                       -                 time-varying
  T_w       = parameters[9]   # Temperature of fusion of ice     K                 273.05              F&MM 1997
  alb_oce   = parameters[10]  # Albedo of ocean                  -                 0.06

  # Get the forcing 
  # ---------------
  forcing_sw = forcing[0]
  forcing_of = forcing[1]

  # Get the state now
  # -----------------
  V_ice_now = state_now[0]
  A_ice_now = state_now[1]
  T_ice_now = state_now[2]
  T_oce_now = state_now[3]

  # Check that the state is physically realistic, complain otherwise
  # ----------------------------------------------------------------
  if V_ice_now < 0.0 or \
     A_ice_now < 0.0 or A_ice_now > 1.0 or  \
     T_ice_now > T_w or T_ice_now <= 0.0 or \
     T_oce_now < T_f or T_oce_now > 100.0 + T_f:
    sys.exit("(sim1D) State unphysical. \n State provided: " + str(state_now) + \
           "\nAllowed values are: \n min: [0.0, 0.0, 0.0+, " + str(T_f) + "] and \n max: [inf, 1.0, " + str(T_w) + ", 373.15]")

  # Exit if the state is not compatible (e.g. ice volume zero but ice concentration not zero)
  # -----------------------------------------------------------------------------------------
  if (np.nanmax((V_ice_now, A_ice_now, T_ice_now)) > 0.0 and np.nanmin((V_ice_now, A_ice_now, T_ice_now)) == 0.0):
    # If at least one of the three variables is zero but at least one is non-zero, then something
    # is wrong
    sys.exit("(sim1D) State incompatibility \n \
                      State provided: " + str(state_now) + "\nIncompatible input concentration, volume and ice surface temperature")

  # Start the physics
  # -----------------
    
  # 1. Update the ocean state and, possibly, create new ice
  # -------------------------------------------------------

  # The net heat supplied to (>0) /lost by (<0) the ocean per unit time is first determined. This consists in a flux directly
  # in the ocean (on an area "1.0 - A_ice_now") and a damped flux through the ice (on an area A_ice_now)
  E_net_oce = (1.0 - A_ice_now) * ((1.0 - alb_oce) * forcing_sw + forcing_of - sigma * T_oce_now ** 4) 
  # Correct for the possible presence of ice
  if V_ice_now > 0.0:
    E_net_oce += A_ice_now      * ((1.0 - alb_ice) * forcing_sw * np.exp(-1.5 * (np.max((V_ice_now / A_ice_now, h_ice_min)) - h_ice_min)))
    # The exponential law follows Fichefet and Morales Maqueda, JGR, 1999

  # That energy is used to cool (<0) or warm (>0) the ocean, depending on the sign.
  # The efficiency at which the ocean changes temperature depends on its specific
  # heat capacity and its mass

  # First, let us determine the amount of energy required to bring the mixed
  # layer to the freezing point, expressed as energy per unit time
  E_f =  (T_f - T_oce_now) * cp_oce * rho_oce * h_oce / dt

  # If the energy balance (E_net_oce) is less than E_f (which is always <=0),
  # then a part is used to create ice
  if E_net_oce < E_f:
    dE = E_f - E_net_oce
    # The seawater temperature is now set to freezing point
    T_oce_new = T_f
    # Sea ice is created in the lead. The volume created per unit time is equal to
    dVdt_ice_accr = dE / Lf_ice
    # That volume is spread on a surface so that the average thickness is at least h_ice_min
    # (If it's more, then the open water is completely filled)
    # For this we first determine the equivalent thicknes obtained if the whole volume is spread
    # on the open water
    h_eq = dVdt_ice_accr / (1.0 - A_ice_now)
    # Then we spread the volume on at least h_ice_min, possibly more
    dAdt_ice_accr = dVdt_ice_accr / np.max((h_ice_min, h_eq))
  else: 
    # Cool/warm water only, without ice production
    dTdt_oce = E_net_oce / (cp_oce * rho_oce * h_oce)
    T_oce_new = T_oce_now + dt * dTdt_oce
    # No ice accretion in that case
    dAdt_ice_accr = 0.0
    dVdt_ice_accr = 0.0

  # 2. Changes in sea ice volume and concentration due to surface melting and bottom accretion/melting
  # --------------------------------------------------------------------------------------------------
  # This step is only necessary if, at the beginning of the time step, ice was present. 
  # The ice possibly created in the last few lines is not taken into account here, as it has already
  # been created during this time step - it does not need to grow even further for now.
  if A_ice_now > 0.0 and V_ice_now > 0.0:
    # Compute conductive heat flux. Here, sea ice thickness is assumed to be uniformly
    # distributed between zero and twice the mean thickness. The
    # conductive flux is therefore boosted with respect to the case of a flat, homogeneous
    # floe. This is because the conductive flux depends on the inverse of thickness.
    # Let us first determine the average thickness of the pre-existing floe
    h_mean_now = V_ice_now / A_ice_now

    # Define an effective conductivity, equal to the actual conductivity boosted to account
    # for the non-flat distribution of thickness. Here, one assumes uniform distribution
    # between 0 and twice the mean thickness.
    # In the absence of snow, the conductive flux is Q = kdT / H where H is the local thickness
    # and kdT is the ice conductivity times the top-bottom temperature difference.

    # The area covered by ice in the thickness interval [h, h + dh]
    # is just dh g(h) = dh * A_ice_now / 2 * h_mean_now (due to the assumption on the ITD).
    # We sum the conductive flux kdT / h over all thickness bins and weight each contribution
    # with the area of ice concerned.
    #
    # For very thin ice (h < h_ice_min), the fluxes go to infinity so we assume a constant
    # flux equal to k dT / h_ice_min (this one covers an area h_min / (2 * h_ice_mean)
    # For ice thicker than h_ice_min, we just integrate dQ = 1 / (2 * h_ice_mean) * k dT / h * dh
    #
    # The resul is k dT / h_ice_mean * (0.5 + 0.5 * ln(2 * h_mean_now / h_ice_min))
    #             = baseline flux * (1.0 + 0.5 * ln(2 * h_mean_now / (e * h_ice_min)))

    boost = 1.0 + 0.5 * np.log(2 * h_mean_now / (np.exp(1) * h_ice_min))
    # The correction factor is by definition > 1.0. For very thin ice 
    # we re-set it to 1.0
    if h_mean_now < np.exp(1) * h_ice_min / 2.0:
      boost = 1.0

    k_ice_eff = k_ice * boost

    # Compute heat imbalance at the top of the floe
    # The net flux on a very thin layer at the top of the ice is:
    #
    # q_net = (1 - alb_ice) * q_sw + q_of - s_bolt * T_ice **4 + q_cond
    #
    # where q_cond = - k_ice_eff * (T_ice - T_f) / h_mean_now (conduction flux)
    #
    # Let's investigate if q_net can equal zero or not. For this we
    # re-write the equation q_net = 0 for the unknown T_ice:
    #
    # 0 = - sigma     * T_ice ** 4 + 
    #       0         * T_ice ** 3 +
    #       0         * T_ice ** 2 + 
    #       k / h_mean_now * T_ice      + 
    #       (1 - alb_ice) * q_sw + q_of - k * T_f / h_mean_now
    # We then solve this equation of the 4th degree using the built-in function
    # of Python.
    A = - sigma
    B = 0.0
    C = 0.0
    D = - k_ice_eff / h_mean_now
    E = (1 - alb_ice) * forcing_sw + forcing_of + k_ice_eff * T_f / h_mean_now
    
    roots = np.roots(np.array((A, B, C, D, E)))

    # We are only interested in the solutions that are not imaginary.
    T_ice_new = np.max(np.real(roots[np.real(roots) == roots]))
  
    # If the new surface temperature (the one that guarantees thermodynamic
    # equilibrium) is larger than the melting point of ice, then the excess
    # of heat must be used to melt the ice at the top, and surface ice temperature
    # must be reset to its melting point.
    if T_ice_new  > T_w:
      # Compute heat imbalance
      q_net_top = - sigma * T_w ** 4 + (- k_ice_eff * (T_w - T_f) / h_mean_now) + (1.0 - alb_ice) * forcing_sw + forcing_of
      # [note!] normally part of the energy must first be used to bring sea ice
      # temperature to freezing point but the heat capacity of sea ice is ignored in this simple model.

      # Compute the rate of surface melting
      dhdt_top = - q_net_top / Lf_ice # Rate of change of thickness at top
     
      # It can happen, due to the numerical resolution of the heat balance equation, 
      # that the growth rate is positive at the top of the ice floe, which is unphysical of course.
      # If it's truly positive then we must stop. Otherwise, if it's slightly positive 
      # (like, if it grows less than a millimeter in one month) we may proceed

      if dhdt_top > 1e-3 / (86400.0 * 30):
        print("(sim1D)     E R R O R : dh/dt|top > 0")
        print("            State before update:")
        print("            h_mean_now = " + str(round(h_mean_now)) + " m")
        print("            T_ice_now = " + str(round(T_ice_now)) + " K")
        print("            Incoming fluxes:")
        print("               Short-wave   ((1 - alb_ice) * q_sw): " + str(round((1 - alb_ice) * forcing_sw)) + " W/m2")
        print("               Other fluxes (                q_of): " + str(round(forcing_of)) + " W/m2")
        print("            Outgoing flux  :")
        print("               Long-wave    (sigma * T ** 4): " + str(round(sigma ** T_ice_now ** 4)) + " W/m2")
        print("            Conduction flux:" + str(- k_ice_eff * (T_w - T_f) / h_mean_now ) + "W/m2")
        sys.exit("(sim1D) E R R O R : dh/dt > 0")            

      # Reset surface temperature to melting point
      T_ice_new = T_w 

    # Else if surface temperature remains below freezing point, nothing happens (ice can't form from the top - snow ice formation is ignored here)
    else:
      dhdt_top = 0.0

    # If new ice was formed, we average the surface temperatures assuming that the newly formed ice is at T_w
    if dAdt_ice_accr > 0.0:
      T_ice_new = (A_ice_now * T_ice_new + dt * dAdt_ice_accr * T_w) / (A_ice_now + dt * dAdt_ice_accr) 

    # 3. Basal accretion or ablation
    # ------------------------------
    # Ice bottom conduction flux based on previous top temperature (one could use
    # the new temperature but I find it more fair to compute all fluxes in terms 
    # of current state variables). Similarly, the flux is here calculated using the
    # old thickness, not the one updated after surface melting.
    q_ice_bot = - k_ice_eff * (T_ice_now - T_f) / h_mean_now 
    # Determine ocean heat flux using bulk parameterization.
    # Follows Goosse and Fichefet, Importance of ice-ocean interactions for the global ocean circulation: A model study, JGR 1999
    # We assume here that the ocean is at rest
    u_ice = 0.01   # Assumed ice velocity, m/s
    stress = rho_oce * 5e-3 * u_ice ** 2 # Ocean-ice stress, with drag coefficient of 0.005 (Goosse and Fichefet, 1999)
    u_star = np.sqrt(stress / rho_oce)   # Friction velocity
    q_oce = rho_oce * cp_oce * 0.006 * u_star * (T_oce_now - T_f)  # The 0.006 is from McPhee, 1992 - see Goosse and Fichefet 1999

    # Net flux at the ice bottom
    q_net_bot = q_oce - q_ice_bot
    # Bottom ice production/ablation rate
    dhdt_bot = - q_net_bot / Lf_ice 

    # CASE 1. ICE IS MELTING
    if (dhdt_bot + dhdt_top) < 0.0:
      # Following Fichefet and Morales Maqueda JGR 1999 we assume that the ice is 
      # a uniform slab except for the update of sea ice concentration, where we 
      # assume that it is uniformly distributed between 0 and twice the mean thickness

      # First let's determine the area of ice ablated. This is the area of ice that 
      # was thinner than - dt * (dhdt_bot + dhdt_top). Since the PDF was 
      # A_ice_now / (2 * h_mean), the rate of ice concentration change due to ablation
      # is:
      dAdt_ice_abl = - A_ice_now * (- (dhdt_bot + dhdt_top)) / (2 * h_mean_now) 

      # The new mean thickness:
      h_mean_new = h_mean_now + dt * (dhdt_bot + dhdt_top)

      # The new volume is computed accordingly, based on the new thickness and the existing area
      dV_ice_topbot = A_ice_now * dt * (dhdt_bot + dhdt_top)
 
    # CASE 2. ICE IS GROWING
    else:
      # Ice has grown. A slab of area A_ice_now times thickness produced is added at the base
      # No lateral accretion is assumed here (this was done above already).
      dAdt_ice_abl = 0.0
      dV_ice_topbot = A_ice_now * dt * (dhdt_bot + dhdt_top)
      h_mean_new = (V_ice_now + dV_ice_topbot) / A_ice_now

  else: # If there was no ice at the beginning of the time step, then only (possibly)
        # ice was created in the leads, but no thermodynamic processes on pre-existing ice occured
    dAdt_ice_abl     = 0.0
    dV_ice_topbot    = 0.0
    h_mean_now       = 0.0
    h_mean_new       = 0.0
    dhdt_top         = 0.0
    dhdt_bot         = 0.0
    
    # If new ice was formed , the surface temperature is set to melting point
    if dAdt_ice_accr > 0.0:
      T_ice_new = T_w
    else:
      T_ice_new = T_ice_now

  # Update the ice volume and concentration
  # ---------------------------------------
  # We update the volume. The ice created by accretion must be there anyway, even if all volume
  # of the thick floe was totally lost
  V_ice_new = np.max((0.0, V_ice_now + dV_ice_topbot)) + dVdt_ice_accr * dt 
  # We do the same for ice concentration
  A_ice_new = np.max((0.0, np.min((A_ice_now + dt * (dAdt_ice_accr + dAdt_ice_abl), 1.0))))
  
  # Set variables to zero/nan if no volume or area 
  # ----------------------------------------------
  if V_ice_new == 0.0 or A_ice_new == 0.0:
    V_ice_new = 0.0
    A_ice_new = 0.0
    T_ice_new = np.nan

  # Debug
  # -----
  #print("")
  #print(":::")
  #print("A_ice_now: " + str(A_ice_now))
  #print("V_ice_now: " + str(V_ice_now) + " m")
  #print("---")
  #print("dA_ice_accr: " + str(dAdt_ice_accr * dt))
  #print("dV_ice_accr: " + str(dVdt_ice_accr * dt) + " m")
  #print("h_mean_now: " + str(h_mean_now) + " m")
  #print("dh_top: " + str(dhdt_top * dt))
  #print("dh_bot: " + str(dhdt_bot * dt))
  #print("h_mean_new: " + str(h_mean_new))
  #print("dA_ice_abl: " + str(dAdt_ice_abl * dt))
  #print("dV_ice_topbot: " + str(dV_ice_topbot) + " m")
  #print("---")
  #print("A_ice_new: " + str(A_ice_new))
  #print("V_ice_new: " + str(V_ice_new) + " m")
  #print(":::")
  #print("")


  # Return the new state
  # --------------------
  state_new = np.empty((4))
  state_new[:] = np.nan
  state_new[0] = V_ice_new
  state_new[1] = A_ice_new 
  state_new[2] = T_ice_new
  state_new[3] = T_oce_new 

  return state_new




# ---------------------------------
# (2) The function to perform a run
# ---------------------------------

def run_sim1D(IC = [1.0, 0.5, -10.0 + 273.15, -1.8 + 273.15], \
              sigma        = 5.67e-8,                         \
              T_f          = -1.8 + 273.15,                   \
              cp_oce       = 4001.896,                        \
              rho_oce      = 1024.458,                        \
              h_oce        = 30.0,                            \
              h_ice_min    = 0.1,                             \
              Lf_ice       = 300.330 * 1e6,                   \
              k_ice        = 2.0344,                          \
              T_w          = 273.05,                          \
              alb_oce      = 0.06,                            \
              alb_ice_mean = 0.774,                           \
              q_sw_mean    = 103.25,                          \
              q_of_mean    = 221.97,                          \
              nt           = 365 * 30,                        \
              dt           = 24.0 * 3600,                     \
              seed         = None,                            \
             ):
  """
  Function run_sim1D
    INPUTS:
      IC = a Numpy array of dimension 4 by 1, providing
           the initial conditions for the sim1D model.
           [volume, area, ice surface temperature, ocean temperature]

      sigma,       |
      T_f,         |
      cp_oce,      |
      rho_oce,     | 
      h_oce,       | ===> SEE THE HELP OF 
      h_ice_min,   | ===> FUNCTION sim1D
      Lf_ice,      |
      k_ice,       |
      T_w,         |
      alb_oce,     |

      alb_ice_mean = annual-mean value for ice albedo
   
      q_sw_mean = annual-mean value for shortwave flux
      q_of_mean = annual-mean value for other fluxes

      nt = number of time steps to integrate forward
      dt = time step in seconds
 
      seed = a seed (integer) to control replicability (the forcing has a random component)

    OUTPUTS:
      A list with two elements:
        - The state as a function of time (Numpy array of size 4 by nt)
        - A datetime object (time axis). By convention the integration is 
          started on January 1st, 1850.
  """

  # PART I. Creation of the forcing 
  # -------------------------------
  day = np.arange(365)

  # Short wave. From Maykut and Untersteiner, 1971; analytical formulation from Notz, PhD thesis, 2004.
  # ---------------------------------------------------------------------------------------------------
  q_sw = 314.0 * np.exp(- 0.5 * ((day + 1 - 164.0) / 47.9) ** 2) + \
         (q_sw_mean - 103.25)

  # Other fluxes. From Maykut and Untersteiner, 1971; analytical formulation from Notz, PhD thesis, 2004.
  # -----------------------------------------------------------------------------------------------------
  q_of = 118.0 * np.exp(- 0.5 * ((day + 1 - 206.0) / 53.1) ** 2) + 179.0 + \
         (q_of_mean - 221.97)

  # Albedo-dependence on day of year. Based on SHEBA campaign. Analytical formulation from Notz, PhD thesis, 2004.
  # --------------------------------------------------------------------------------------------------------------
  albedo_ice = 0.914 - 0.431 / (1.0 + ((day + 1.0 - 207) / 44.5) ** 2) + \
         (alb_ice_mean - 0.774)
               

  # Generate perturbation to add to the mean fluxes (interannual variability)
  # -------------------------------------------------------------------------
  if seed is not None:
    np.random.seed(seed)

  # We generate auto-correlated noise added on top of interannual variability.
  q_pert_of = np.empty((nt))
  q_pert_of[:] = np.nan
  q_pert_of[0] = 3.0 * np.random.randn()

  q_pert_sw = np.empty((nt))
  q_pert_sw[:] = np.nan
  q_pert_sw[0] = 4.0 * np.random.randn()

  for jt in range(1, nt):
    q_pert_of[jt] = 0.90 * q_pert_of[jt - 1] + 3.0 * np.random.randn()
    q_pert_sw[jt] = 0.50 * q_pert_sw[jt - 1] + 4.0 * np.random.randn()


  # PART II. Integration
  # --------------------
  state = np.empty((4, nt))
  state[:] = np.nan
  state[:, 0] = IC

  for jt in range(1, nt):
    # Determine day of year.
    doy = np.int(np.mod(jt * dt, 365.0 * 24.0 * 3600.0) / (24.0 * 3600.0))

    # Determine the forcing for that day
    forcing = [ np.max((q_sw[doy]  + q_pert_sw[jt], 0)),       \
                np.max((q_of[doy]  + q_pert_of[jt], 0)),       \
              ]

    parameters = [sigma, \
                  T_f,   \
                  cp_oce, \
                  rho_oce, \
                  h_oce, \
                  h_ice_min, \
                  Lf_ice, \
                  k_ice, \
                  albedo_ice[doy], \
                  T_w,         \
                  alb_oce]

    state[:, jt] = sim1D(state[:, jt - 1], parameters, forcing, dt = dt)


  base = datetime.datetime(1850, 1, 1)
  t = list()
  jt = 0
  jtt = 0
  while jt < nt:
    tmp = base + datetime.timedelta(days = jtt * dt / (24.0 * 3600) ) 
    if tmp.month != 2 or tmp.day != 29:
      t.append(tmp)
      jt += 1
      jtt += 1
    else:
      jtt += 1

  return [state, t]



if __name__ == '__main__':

  [state, t] = run_sim1D(seed = 1)
  
  offset = [0.0, 0.0, -273.15, -273.15]
  scalef = [1.0, 100.0, 1.0, 1.0]
  lab = ["Volume [m]", "Concentration [%]", "Ice temperature [C]", "Mixed layer\ntemperature [C]"]
  plt.figure(figsize = (8, 12))
  for j in range(4):
    plt.subplot(4, 1, j + 1)
    plt.plot(t, scalef[j] * state[j, :] + offset[j], lw = 2)
    plt.xlim(t[0], t[-1])
    plt.grid()
    plt.ylabel(lab[j])

  plt.tight_layout()
  plt.savefig("./figs/tmp.png", dpi = 300)
  print("Figure ./figs/tmp.png printed")
  plt.clf()  

