#!/usr/bin/env python
#
#
# Copyright (c) 2013, Unbounded Robotics Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

from math import sin, cos, pi, sqrt


####################################################################################################
#
# Implements Field-oriented control (FOC) for sinusoidally would permanent magnet brushless motor.
# Field-oriented control is also refered to as Vector Control
# 
# There a few reference descibing FOC for sinusoidally wound brushless permentent magnet motors.
# Many overly complicate the issues, or focus on new research area.  
# The best reference for the FOC algorithm used in this code is a app-note from Texas Instruments
#
# Sensorless Field Oriented Control of Multiple Permanent Magnet Motors
#  Bilal Akin
#  Manish Bhardwaj
#
# Fig 7 in the app-note shows a basic FOC control scheme.
#
#
# FOC block diagram
#  
#           +     +------------+ Vd_ref   +-----------+            +--------+        +---------+
# Id_ref --->O--> | PI         | -------> |           | Valpha_ref |        | DutyA  |         |
#           -^    | Controller |          |           | ---------> |        | -----> |         |
#            |    +------------+          |           |            | Quasi  |        | 3-phase |
#            |                            | Inverse   |            | Space- | DutyB  | Inverter|
#    Id -----'                            | Park      | Vbeta_ref  | Vector | -----> |         |
#                                         | Transform | ---------> | PWM    |        |         |
#           +     +------------+ Vq_ref   |           |            |        | DutyC  |         |
# Iq_ref --->O--> | PI         | -------> |           |            |        | -----> |         |
#           -^    | Controller |          +-----------+            +--------+        +---------+
#            |    +------------+                ^                                       | | |
#            |                                  |                                       | | |
#    Iq -----'                                  +-- rotor angle +                       | | |
#                                                     rotor velocity * delay            | | |
#                                                                                       | | |
#                                                                                       | | |
#                                              +--- rotor angle                         | | |
#                                              |                                        | | |
#                                              v                                        | | |
#                                        +-----------+          +-----------+    I_A    | | |
#                                   Id   |           | I_alpha  |           | <---------( | | 
#                                <------ |  Park     | <------- | Clark     |    I_B    | | |
#                                   Iq   | Transform | I_beta   | Transform | <-----------( |
#                                <------ |           | <------- |           |    I_C    | | |
#                                        |           |          |           | <-------------(
#                                        +-----------+          +-----------+           |_| | 
#                                                                                      /     \
#                                                                                     | Motor |
#                                                                                      \_____/
#                                                                                        | |
#                                                                        rotor angle <---+ | 
#                                                                                          | 
#                                                                   rotor velocity <-------+ 
#
####################################################################################################


####################################################################################################
# Torque Constant for Sinusoidal Commutation:
#    Most motor torque specs are based on block commutation even though the motor is 
#    sinusoidally wound. The question is how to modify the block-commutation torque constant 
#    and torque limits for use with field-oriented control.
#    Luckily, maxon has hard to find app-note about how to map parameters based 
#    on block commutation to parameters used for sinusoidal commutation:
#      "Parameters of maxon-EC-Motors Driven By Servo Amplifiers With Sinus-Commutation"
#
# The app note comes up with the equation (16) :
#   To = pi/(2*sqrt(3)) * Kb * I_s
#
# In the equation :
#   To is the output torque
#   Kb is the torque constant given on the motor datasheet that was determined using block communtation of the motor
#   I_s the the magnitude of the current 
#
# However it is somewhat unclear how I_s maps to the Iq value that we use for torque. 
# Based on running the rest the the equations in the app-note:
#   iq = sqrt(3/2)*I_s
#     or
#   I_s = sqrt(2/3)*I_q
#
#
# This gives us the final formula
#   To = pi/(2*sqrt(3)) * Kb * I_s
#   To = pi/(2*sqrt(3)) * Kb * sqrt(2/3)*I_q
#   To = pi*sqrt(2)/(2*sqrt(3)*sqrt(3)) * Kb * Iq
#   To = pi*sqrt(2)/(2*3) * Kb * Iq
#   To = pi/(sqrt(2)*3) * Kb * Iq
#   To = 0.74 * Kb * iq
#
#
# As far the the max torque torque of the motor, three equations are used
#   To = pi/(2*sqrt(3)) * Kb * I_s
#   Is_max = 2/sqrt(3) * Ib_max
#   To_max(b) = Ib_max * Kb   or (Ib_max = To_max(b) / Kb)
#
#  Putting these two equations together
#   To_max(s) = pi/(2*sqrt(3)) * Kb * Is_max
#   To_max(s) = pi/(2*sqrt(3)) * Kb * 2/sqrt(3) * Ib_max
#   To_max(s) = pi/(2*sqrt(3)) * Kb * 2/sqrt(3) * To_max(b) / Kb
#   To_max(s) = pi/(2*sqrt(3)) * 2/sqrt(3) * To_max(b)
#   To_max(s) = pi/(2*sqrt(3)) * 2/sqrt(3) * To_max(b)
#   To_max(s) = pi/3 * To_max(b)
#   To_max(s) = 1.047 * To_max(b)
#
#  Conclusion :
#     Torque = pi/(sqrt(2)*3) * Kb * Iq   where Kb is value provided by Maxon for block commutation
#     MaxTorque = pi/3 * MaxTorque(b)     where MaxTorque(b) is max toruqe value provided by Maxon 
####################################################################################################


####################################################################################################
# Space Vector Modulation 
#   SVM generates output voltage by treating both output voltage and possible 
#   switch combinations as vectors.  
#   The output voltage is generated by modulating between:
#     * The two switch vectors that are adjacent to the output vector
#     * One or both of the null vectors (all switches high or all switches low)
#  
#   Advantages of SVM 
#     The problem with treating all output voltages separately, is that it limit the possible 
#     voltage difference between the output phases.  This limits that max motor speed.
#     SVM treats 3 PWMs as a single unit, allow a 15.5% higher output voltage.
#
#   Choice of Null Vector
#     With SVM there is open choice of which null vector to use and when.  
#     In our case, it seems like the best to always use the V0 null vector (all switches output low)
#     The advantages of doing this is that one of the phases never switches (is always low)
#     during a complete PWM cycle.  This reduces the switch power loss, EMI, etc...
#     The same reduction in switching is also possible by always choosing the 
#     V7 null vector (all switch outputs high).  However this is a bad idea since the phase 
#     current can only sampled when the phase is low (on our design).
#
#   Why not SVM
#     The book implementation uses (or implies using) arctan() to determine the angle of the 
#     desired output phase.  Then it uses this angle to find a appropriate entry in a table 
#     of formulas.  The formulas determine duty for each PWM and are pretty simple.
#     The TI app-note suggests a method of indexing into the table of formulas without using 
#     arctan().  This method basically performs a (modified) inverse Clark transform and uses
#     some trick with sign() to come up with a table index.
#    
#   A better SVM implementation (for always Null=V0)
#       1. First use normal inverse Clark transform to determine amplitude of the 3 output voltages
#       2. Determine minimum of value of three output voltages
#       3. Subtract minimum value from all 3 output voltages 
#       4. ?
#       5. Profit!
#     This should use fewer operations than even the TI implementation. 
#     It is also probably less error prone to implement (no table of formulas).
#     However, it can not produce the other SVM variations.
#
####################################################################################################

####################################################################################################
# Phase Alignment
#   To properly drive a 3-phase brushless motor we need a acurrate rotor position measurement.
#   When an absolute encoder is used it is easy to get a good position measurement.
#   However, the absolute encoder must be align to the phases of the motor at least once before 
#   it is used. 
#
#   When encoder is *aligned* to rotor, a encoder angle of 0 should correpond to rotor direct axis
#   being align with motor phase A.
#
#   A very naive method of aligning motor would be to drive a strong field in the direction of 
#   phase A.  Since the torque = K*I*sin(stator_current_angle - rotor_angle).  
#   The torque will be 0 with rotor_angle == stator_current_angle.   
#   However, the torque is also zero when (rotor_angle - stator_current_angle) = 180 degrees.
#   Which means the rotor is at an unstable equalibrium point. 
#   Also in the presense of static friction or external load, the motor might not make it 
#   all the way to the zero position.
#
#   Another method of alignment would be to drive the motor externally with allowing current 
#   to flow on motor phase.  This should produce sine back-emf voltage on all three phases.
#   The zero-crossing of the back-emf the voltage on phase A is also the direct axis of the 
#   motor.  By determining where the this zero-crossing occurs, the angle of the rotor direct axis
#   can be found.  
#   Unfortunately this method can only be performed on a test-stand because an external torque 
#   is needed to drive the motor.  The motor cannot drive itself, the current used would 
#   change the phase voltage and skew the zero-crossing position because of the phase resistance
#   and inductance.
#
#   For an motor that must drive itself, the zero-crossing can be found by subtracting out the 
#   voltage contributions of the current by using know values of motor resistance and
#   inductance.  If a rough location of the rotor direct axis is known, it is even possible 
#   to iteratively find the direct axis while the motor is moved at a high speed.
#   This is explained below.
#   
# Iterative alignement
#   Since the FOC control algorithm already breaks motor state into direct and quadrature values,
#   is easiest to work with these values instead of sinusoidal values.
#   The "Vector Control of Three-Phase AC Machines" provides the following voltage equation 
#   for the permanant magent brushless motor in Q-D coordinates is given in section 3.3.1
#
#   Equation 3.60;
#     Vd = Rs * Id + Ld * dId/dt - el_velocity*Lq*Iq
#     Vq = Rs * Iq + Lq * dIq/dt + el_velocity*Ld*Id + el_velocity*pole_flux 
#
#   Assume that the actual rotor angle has an error A.  Then actual values of 
#   of Vd and Vq the controller see will be slighly wrong.  
#   To determine values of Vd and Vq the controller calculates (Vd' and Vq'), 
#   we would rotate the vectors by A using a rotation matrix:
#     Vd' = Vd * cos(A) - Vq * sin(A)
#     Vq' = Vd * sin(A) + Vq * cos(A)
# 
#   The controllers idea of Id and Iq is also wrong:
#     Id' = Id * cos(A) - Iq * sin(A)
#     Iq' = Id * sin(A) + Iq * cos(A)
#  
#   Now if we assume that A is small then cos(A) ~= 1 and sin(A) ~= A 
#     Vd' = Vd - Vq * A
#     Vq' = Vd * A + Vq 
#     Id' = Id - Iq * A
#     Iq' = Id * A + Iq
#
#   If we write out equation for Vd'
#     Vd' =   [Rs * Id + Ld * dId/dt - el_velocity*Lq*Iq] - A*[Rs * Iq + Lq * dIq/dt + el_velocity*Ld*Id + el_velocity*pole_flux]
#     Vq' = A*[Rs * Id + Ld * dId/dt - el_velocity*Lq*Iq] +   [Rs * Iq + Lq * dIq/dt + el_velocity*Ld*Id + el_velocity*pole_flux]
#
#   Here's what controller's estimate of Vd
#     Vd_est = Rs * Id' + Ld * dId'/dt - el_velocity*Lq*Iq'
#
# 
#   If we subtract controller estimate of Vd from measured Value of Vd we get an error
#     error = Vd' - Vd_est 
#     error = [Rs * Id + Ld * dId/dt - el_velocity*Lq*Iq] - A*[Rs * Iq + Lq * dIq/dt + el_velocity*Ld*Id + el_velocity*pole_flux] - [Rs * Id' + Ld * dId'/dt - el_velocity*Lq*Iq']
#
#   Now replace Id and Iq with Id' and Iq'
#     Id = (Id' + A*Iq')
#     Iq = (-A*Id' + Iq')
#
#     error = [Rs *       Id      + Ld *      dId/dt       - el_velocity*Lq*Iq]             - A*[Rs *        Iq      + Lq * dIq/dt             + el_velocity*Ld*Id            + el_velocity*pole_flux] - [Rs * Id' + Ld * dId'/dt - el_velocity*Lq*Iq']
#     error = [Rs * (Id' + A*Iq') + Ld * d(Id' + A*Iq')/dt - el_velocity*Lq*(-A*Id' + Iq')] - A*[Rs * (-A*Id' + Iq') + Lq * d(-A*Id' + Iq')/dt + el_velocity*Ld*(Id' + A*Iq') + el_velocity*pole_flux] - [Rs * Id' + Ld * dId'/dt - el_velocity*Lq*Iq']
#     error = Q1 + Q2 + Q3
#
#     Q1 = Rs * (Id' + A*Iq') + Ld * d(Id' + A*Iq')/dt - el_velocity*Lq*(-A*Id' + Iq')
#     Q1 = Rs*Id' + A*Rs*Iq' + Ld*dId'/dt + A*Ld*Iq'/dt + A*el_velocity*Lq*Id' - el_velocity*Lq*Iq'
#
#     Q2 = - A* [Rs * (-A*Id' + Iq') + Lq * d(-A*Id' + Iq')/dt + el_velocity*Ld*(Id' + A*Iq') + el_velocity*pole_flux]
#     Q2 = - A* [-A*Rs*Id' + Rs*Iq' - A*Lq*dId'/dt + Lq*dIq'/dt + el_velocity*Ld*Id' + A*el_velocity*Ld*Iq' + el_velocity*pole_flux]
#     Q2 = A*A*Rs*Id' - A*Rs*Iq' + A*A*Lq*dId'/dt - A*Lq*dIq'/dt - A*el_velocity*Ld*Id' - A*A*el_velocity*Ld*Iq' - A*el_velocity*pole_flux
#
#     Q3 = - [Rs * Id' + Ld * dId'/dt - el_velocity*Lq*Iq']
#     Q3 = -Rs*Id' - Ld*dId'/dt + el_velocity*Lq*Iq'
#
#     Q1 + Q3 = Rs*Id' + A*Rs*Iq' + Ld*dId'/dt + A*Ld*Iq'/dt + A*el_velocity*Lq*Id' - el_velocity*Lq*Iq' - Rs*Id' - Ld*dId'/dt + el_velocity*Lq*Iq'
#     Q1 + Q3 = A*Rs*Iq' + A*Ld*Iq'/dt + A*el_velocity*Lq*Id'
#     Q1 + Q3 = A*Rs*Iq' + A*Ld*Iq'/dt + A*el_velocity*Lq*Id'
#
#     error = Q1+Q2+Q3 = A*Rs*Iq' + A*Ld*Iq'/dt + A*el_velocity*Lq*Id' + A*A*Rs*Id' - A*Rs*Iq' + A*A*Lq*dId'/dt - A*Lq*dIq'/dt - A*el_velocity*Ld*Id' - A*A*el_velocity*Ld*Iq' - A*el_velocity*pole_flux
#
#   Now assume Lq == Ld == L:
#     error = A*Rs*Iq' + A*Ld*Iq'/dt + A*el_velocity*Lq*Id' + A*A*Rs*Id' - A*Rs*Iq' + A*A*Lq*dId'/dt - A*Lq*dIq'/dt - A*el_velocity*Ld*Id' - A*A*el_velocity*Ld*Iq' - A*el_velocity*pole_flux
#     error = A*Rs*Iq' + A*L*Iq'/dt + A*el_velocity*L*Id' + A*A*Rs*Id' - A*Rs*Iq' + A*A*L*dId'/dt - A*L*dIq'/dt - A*el_velocity*L*Id' - A*A*el_velocity*L*Iq' - A*el_velocity*pole_flux
#     
#   Pull out A
#     error = A * [Rs*Iq' + L*Iq'/dt + el_velocity*L*Id' + A*Rs*Id' - Rs*Iq' + A*L*dId'/dt - L*dIq'/dt - el_velocity*L*Id' - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#
#   Simplify
#     error = A * [Rs*Iq' + L*Iq'/dt + el_velocity*L*Id' + A*Rs*Id' - Rs*Iq' + A*L*dId'/dt - L*dIq'/dt - el_velocity*L*Id' - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#     error = A * [       + L*Iq'/dt + el_velocity*L*Id' + A*Rs*Id'          + A*L*dId'/dt - L*dIq'/dt - el_velocity*L*Id' - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#     error = A * [       + L*Iq'/dt                  + A*Rs*Id'          + A*L*dId'/dt - L*dIq'/dt                  - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#     error = A * [                                   + A*Rs*Id'          + A*L*dId'/dt                              - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#     error = A * [ A*Rs*Id' + A*L*dId'/dt - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#
#   Now assume controller is doing decent of keeping both Id' and dId' near 0:
#     error = A * [ - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#
#   The controller can calculate the error value, but we want to determine A  
#     A = error / [ - A*el_velocity*L*Iq' - el_velocity*pole_flux]
#     A = -error / [  A*el_velocity*L*Iq' + el_velocity*pole_flux]
#     A = -error / [el_velocity*(A*L*Iq' + pole_flux)]
#
#   Unfortunately, there is an A-term on the right side, so formula will not calculate the correct value of A.  
#   We could wait for Iq term to happen to be zero before calculating A.  If Iq is 0 then the A can be calculated
#   exactly with formula:
#     A = -error / (el_velocity*pole_flux)
#
#   However, Iq is not always guarenteed to be zero.  
#   Even if Iq was 0, it is measured value will non-zero noise and offset.
#   
#   If instead of calculating A in a signle step, the value if estimated value of (A') is use to offset
#   the angle of the motor the next cycle:
#      First calculate estimate of A (A'):
#         A' = -error / (el_velocity*pole_flux)
#      Offset angle for next controller cycle by estimate of A.  This basically changes A 
#         A_next = A - A'
#      If the estimated value (A') is close to real value of A, then abs(A_next) < abs(A) and 
#      the value of A will converge to 0.
# 
#   The about iterative algorithm should work, unless abs(A-A') > abs(A).
#   Writen another way, the algorihtm should work if  abs((A-A')/A) < 1
#      Actual value of A:
#        A = -error / [el_velocity*(A*L*Iq' + pole_flux)]
#        A = -error/el_velocity / (A*L*Iq' + pole_flux)
#      Estimate of of
#        A' = -error/el_velocity / pole_flux
#      Difference between both values (A-A')
#        A = -error / [el_velocity*(A*L*Iq' + pole_flux)] - -error/el_velocity / pole_flux
#        A-A' = -error/el_velocity * [ 1/(A*L*Iq' + pole_flux) - 1/pole_flux ]
#      A-A'/A
#        A-A'/A = [ 1/(A*L*Iq' + pole_flux) - 1/pole_flux ] * (A*L*Iq' + pole_flux)
#        A-A'/A = [ 1 - (A*L*Iq' + pole_flux)/pole_flux ]
#        A-A'/A = [ 1 - (A*L*Iq' + pole_flux)/pole_flux ]
#        A-A'/A = [ 1 - (A*L*Iq' + pole_flux)/pole_flux ]
#        A-A'/A = [ 1 - (A*L*Iq'/pole_flux + 1) ]
#        A-A'/A = [ 1 - (A*L*Iq'/pole_flux + 1) ]
#        A-A'/A = A*L*Iq'/pole_flux
#
#   The alogrithm with converge when abs(A*L*Iq'/pole_flux) is < 1
#   Since the worst case A is 180 degrees = pi radians ~= 3. The 
#   The ratio: 
#      abs(L*Iq'/pole_flux) < 1/3.
#   or 
#      abs(pole_flex/(L*Iq)) > 3
#
# Implementation:
#   First find error between Vd and what it should be 
#     Vd_error = Vd - (Rs * Id + Ld * dId/dt - el_velocity*Lq*Iq)
#  
#   To simple assume both Id and dId/dt is small.
#   Note: dId/dt is not small if Iq (or command Iq) is changing often
#     Vd_error = Vd + el_velocity*Lq*Iq
#
#   Estimate backemf voltage 
#     backemf_voltage = (el_velocity*backemf_constant)
#
#   Estimate offset error (pole_flux = backemf_constant * pole_pairs)
#     est_offset_error = -Vd_error / (el_velocity*backemf_constant*pole_pairs)
#     est_offset_error = -Vd_error / (el_velocity*pole_pairs*backemf_constant)
#     est_offset_error = -Vd_error / (velocity*backemf_constant)
#     est_offset_error = -Vd_error / backemf_voltage
#
#   If backemf voltage if above a certain threshold, update offset_error with estimate
#     offset_error += alpha * est_offset_error
#    




def alignmentRatio(motor):
    """ Calculates ratio of motor (backemf*poles_pole_pair) to (L*Iq_max)
    This ratio determines how well interative alignment routie will converge.
    If ratio is high (> 100) the value will converge easily.  
    If ratio is 3 or less, then algorithm might not converge when both
    intial angle offset is large and commanded torque is high
    """
    # using max torque determine maximum value of Iq
    
    # torque_constant = motor.torque_constant * sqrt(2)*3/pi
    # max_iq = motor.max_torque * pi/3 / torque_constant
    # L = motor.terminal_inductance / sqrt(2)    
    # max_iq * L = motor.max_torque * pi/3 / motor.torque_constant * sqrt(2)*3/pi * motor.terminal_inductance / sqrt(2)
    # max_iq * L = motor.max_torque * motor_terminal_inductance / motor.torque_constant
    #max_iq = motor.max_torque * motor_terminal_inductance / motor.torque_constant
    
    # backemf_constant = motor.torque_constant # Nm/A = V/(rad/s)
    # (backemf_constant*motor.pole_pairs) / (L*max_iq) =
    # (motor.torque_constant*motor.pole_pairs) / (motor.max_torque * motor_terminal_inductance / motor.torque_constant)
    # (motor.torque_constant**2 * motor.pole_pairs) / (motor.max_torque * motor.terminal_inductance)

    ratio = (motor.pole_pairs * motor.torque_constant**2) / (motor.max_torque * motor.terminal_inductance)
    print "Motor alignment ratio is ", ratio
    return ratio    




class PI_Controller:
    """ Specialized PI controller for FOC control, limits both output 
    and I-term to range of -1 to 1
    """
    def __init__(self, gains):
        p_gain, i_gain = gains
        self.p_gain = p_gain
        self.i_gain = i_gain
        self.t = [0.0]
        self.i_term = [0.0]
        self.output = [0.0]

    def update(self, t, dt, error):
        i_term = self.i_term[-1]
        i_term += error*self.i_gain*dt
        if i_term > 1.0: 
            i_term = 1.0        
        elif i_term < -1.0: 
            i_term = -1.0
        output = error*self.p_gain + i_term
        if output > 1.0:
            output = 1.0
        elif output < -1.0:
            output = -1.0
        self.output.append(output)
        self.i_term.append(i_term)
        self.t.append(t)
        return output


class FOC_Controller:
    def __init__(self, dt, motor, d_gains, q_gains, supply_current_limit):
	self.name = "Python FOC"
        self.dt = dt
        self.motor = motor
        self.supply_current_limit = supply_current_limit

        self.d_ctrl = PI_Controller(d_gains)
        self.q_ctrl = PI_Controller(q_gains)
        
        self.t = [0.0]
        self.measured_torque = [0.0]
        self.measured_id = [0.0]
        self.measured_iq = [0.0]

        self.input_target_torque = [0.0] # what torque was requested from higher level
        self.target_torque = [0.0] #what torque was actually limited to
        self.target_id = [0.0]
        self.target_iq = [0.0]
        
        self.cmd_a = [0.0]
        self.cmd_b = [0.0]
        self.cmd_c = [0.0]

        self.cmd_d = [0.0]
        self.cmd_q = [0.0]
        self.cmd_q_filt = [0.0]
        self.cmd_vd = [0.0]
        self.est_vd = [0.0]
        self.error_vd = [0.0]
        self.backemf_voltage = [0.0]
        self.est_el_angle_offset = [0.0]
        self.el_angle_offset = [0.0]

        self.torque_limit = [0.0]

    def reset(self, velocity, supply_voltage):
        self.d_ctrl.i_term[-1] = 0.0        
        self.q_ctrl.i_term[-1] = -velocity * self.motor.torque_constant / supply_voltage * pi/3.0
        self.cmd_q_filt[-1] = 0.0

    def update(self, t, target_torque, position, velocity, ia, ib, ic, supply_voltage):
        dt = self.dt 
        motor = self.motor
        self.input_target_torque.append(target_torque)

        # Conversion from torque to quadrature current
        #     Torque =  pi/(sqrt(2)*3) * Kb * Iq
        #     Iq = Torque/Kb * sqrt(2)*3/pi        
        torque_constant = self.motor.torque_constant * 0.7404804896930609 # sqrt(2)*3/pi

        # Have dynamic limit on torque based on supply current limit (aka a power limit)
        torque_limit = 0.0
        if self.supply_current_limit > 0.0:
            cmd_q_filt = self.cmd_q_filt[-1]
            if abs(cmd_q_filt) > 1e-3:
                sqrt2 = 1.4142135623730951 # sqrt(2)
                torque_limit = -torque_constant*self.supply_current_limit*sqrt2/cmd_q_filt
                if torque_limit > 0.0:
                    target_torque = min(target_torque, torque_limit)
                else:
                    target_torque = max(target_torque, torque_limit)

        # Always enforce motor torque limit for robot
        if target_torque > motor.max_torque:
            target_torque = motor.max_torque
        elif target_torque < -motor.max_torque:
            target_torque = -motor.max_torque

        # Clarke Transform
        #  map 3-phase measured currents (ia,ib,ic) into two phase currents i_alpha, i_beta
        #  i_alpha is aligned with ia axis
        #  i_beta is 90 degrees from ia axis
        #  The sqrt(2/3) scaling factor keeps the magnitude of the vectors equal
        #   ia**2 + ib**2 + ic**2 == i_alpha**2 + i_beta**2
        # Constants (for later)
        #
        s30 = 0.5  # sin(30)
        c30 = 0.8660254037844387 # cos(30)
        sqrt2d3 = 0.816496580927726 # sqrt(2/3)
        i_alpha = sqrt2d3 * (ia  - s30*ib - s30*ic)
        i_beta  = sqrt2d3 * (      c30*ib - c30*ic)

        # Motor electrial angle, 
        #  the number of electical cycles for every mechanical rotor cycle depends 
        #  on the number of motor poles-pairs.
        pole_pairs = self.motor.pole_pairs
        el_angle = position * pole_pairs - self.el_angle_offset[-1]
	

        # Park Tranform
        #  use rotor electrical angle to map i_alpha, i_beta currents into 
        #  rotor-oriented values (i_d, i_q).
        #  i_d = direct axis current
        #  i_q = quadrature axis current
        #  this is basicaly a rotatation matrix multiplication where by -el_angle
        s1 = sin(el_angle)
        c1 = cos(el_angle)
        measured_id =   c1*i_alpha + s1*i_beta
        measured_iq =  -s1*i_alpha + c1*i_beta
    
        # using iq to determine actual torque production 
        # assume id is small enough that is does not effect magnetic field
        measured_torque = -measured_iq * torque_constant

        # Determine target current values
        #   direct current : should always be 0 since we don't want any field weakening
        #   quadrature current : should be proportional to the motor torque
        target_id = 0.0
        target_iq = -target_torque / torque_constant

        # Run PI controllers
        cmd_d = self.d_ctrl.update(t, dt, target_id - measured_id)
        cmd_q = self.q_ctrl.update(t, dt, target_iq - measured_iq)

        # perform circle limit on d/q commands
        # this limits the magnitude of the d/q command to 1.0 while keeping the angle the same
        mag = sqrt(cmd_d*cmd_d + cmd_q*cmd_q)
        if mag > 1.0:
            cmd_d /= mag
            cmd_q /= mag

        # Keep a filtered cmd_q value for supply current limiting
        cmd_q_filt = self.cmd_q_filt[-1]
        cmd_q_filt += (cmd_q - cmd_q_filt) * 0.25

        # Project motor angle into future using time delay and current velocity
        #   delay is based time between measurement of  current values to execution of new commands
        #   this is about 1-2 PWM cycles
        #   this is not absolutely needed but improves performace at higher speeds
        el_velocity = velocity * self.motor.pole_pairs
        el_angle2 = el_angle + (el_velocity * dt) * 1.32

        # Map d/q commands back to stator referenced values (inverse park transform)
        s2 = sin(el_angle2)
        c2 = cos(el_angle2)
        cmd_alpha =  c2*cmd_d - s2*cmd_q
        cmd_beta  =  s2*cmd_d + c2*cmd_q
        
        # Map alpha/beta commands to 3-phase commands (inverse Clarke transform)
        inv_sqrt3 = 0.5773502691896258 # 1/sqrt(3)
        cmd_a = inv_sqrt3 * (      cmd_alpha                )
        cmd_b = inv_sqrt3 * ( -s30*cmd_alpha + c30*cmd_beta )
        cmd_c = inv_sqrt3 * ( -s30*cmd_alpha - c30*cmd_beta )

        # Now perform quasi- SVM (space vector modulation)
        #   See notes about SVM
        cmd_min = min(min(cmd_a, cmd_b), cmd_c)
        cmd_a -= cmd_min
        cmd_b -= cmd_min
        cmd_c -= cmd_min

        # Determine the difference between the commanded value of direct voltage 
        # and what is should be based on velocity and current
        cmd_vd = supply_voltage*self.cmd_d[-1] 
        L = self.motor.terminal_inductance / sqrt(2)
        est_vd = - el_velocity * measured_iq * L 
        error_vd = cmd_vd - est_vd

        # Update estimate of el_angle_offset
        backemf_constant = self.motor.torque_constant
        backemf_voltage = velocity * backemf_constant
        if abs(backemf_voltage) > 1.0:
            est_el_angle_offset = -error_vd / backemf_voltage
            alpha = 0.01
            el_angle_offset = self.el_angle_offset[-1] + alpha * est_el_angle_offset 
        else:
            est_el_angle_offset = 0.0
            el_angle_offset = self.el_angle_offset[-1]

        # record values
        self.t.append(t)
        self.measured_torque.append(measured_torque)
        self.measured_id.append(measured_id)
        self.measured_iq.append(measured_iq)
        self.target_torque.append(target_torque)
        self.target_id.append(target_id)
        self.target_iq.append(target_iq)
        self.cmd_a.append(cmd_a)
        self.cmd_b.append(cmd_b)
        self.cmd_c.append(cmd_c)

        self.cmd_d.append(cmd_d)
        self.cmd_q.append(cmd_q)
        self.cmd_q_filt.append(cmd_q_filt)
        self.cmd_vd.append(cmd_vd)
        self.est_vd.append(est_vd)
        self.error_vd.append(error_vd)
        self.backemf_voltage.append(backemf_voltage)
        self.est_el_angle_offset.append(est_el_angle_offset)
        self.el_angle_offset.append(el_angle_offset)

        self.torque_limit.append(torque_limit)

        return (cmd_a, cmd_b, cmd_c)




class FOC_CWrapper:
    """ Wrapper around C++ implementation of FOC-algorithm 
    Makes interface to C++ functions almost exactly the same as the python implemenation.
    """
    def __init__(self, dt, motor, d_gains, q_gains, supply_current_limit):
	self.name = "C++ FOC"
        # I know this is frowned upon, but we don't want to import foc_py
        # if we are not using it -- especially since it needs to be compiled first
        import foc_py 
        foc = foc_py.FOC()
        self.foc = foc
        foc.enableIncrementalAlignment(True)

        # adjust I gains by dt because C-implementation of FOC d/n do this        
        foc.setDGains(d_gains[0], d_gains[1]*dt) 
        foc.setQGains(q_gains[0], q_gains[1]*dt)
        foc.torque_constant = motor.torque_constant
        foc.max_torque = motor.max_torque
        foc.terminal_inductance = motor.terminal_inductance
        foc.offset_alignment_alpha = 0.01
        foc.min_backemf_for_offset_alignment = 1.0
        foc.pole_pairs = motor.pole_pairs
        foc.supply_current_limit = supply_current_limit

        foc.foc_timestep = dt
        foc.pole_pairs = motor.pole_pairs

        self.cmds_v3 = foc_py.FOC_Vector3()
        
        self.t = [0.0]
        self.input_target_torque = [0.0]
        self.measured_torque = [0.0]
        self.measured_id = [0.0]
        self.measured_iq = [0.0]

        self.target_torque = [0.0]
        self.target_id = [0.0]
        self.target_iq = [0.0]
        
        self.cmd_a = [0.0]
        self.cmd_b = [0.0]
        self.cmd_c = [0.0]

        self.cmd_d = [0.0]
        self.cmd_q = [0.0]
        self.cmd_vd = [0.0]
        self.est_vd = [0.0]
        self.error_vd = [0.0]
        self.backemf_voltage = [0.0]
        self.est_el_angle_offset = [0.0]
        self.el_angle_offset = [0.0]
        self.torque_limit = [0.0]

    def reset(self, velocity, supply_voltage):
        self.foc.reset(velocity, supply_voltage)

    def update(self, t, target_torque, position, velocity, ia, ib, ic, supply_voltage):
        foc = self.foc

        foc.update(self.cmds_v3, target_torque, position, velocity, ia, ib, ic, supply_voltage)
        cmd_a = self.cmds_v3.a
        cmd_b = self.cmds_v3.b
        cmd_c = self.cmds_v3.c


        self.t.append(t)
        self.input_target_torque.append(target_torque)
        self.measured_torque.append(foc.measured_torque)
        self.measured_id.append(foc.measured_id)
        self.measured_iq.append(foc.measured_iq)
        self.target_torque.append(target_torque)
        self.target_id.append(foc.target_id)
        self.target_iq.append(foc.target_iq)
        self.cmd_a.append(cmd_a)
        self.cmd_b.append(cmd_b)
        self.cmd_c.append(cmd_c)

        self.cmd_d.append(foc.cmd_d)
        self.cmd_q.append(foc.cmd_q)
        self.cmd_vd.append(foc.cmd_vd)
        self.est_vd.append(foc.est_vd)
        self.error_vd.append(foc.error_vd)
        self.backemf_voltage.append(foc.backemf_voltage)
        self.est_el_angle_offset.append(foc.est_el_angle_offset)
        self.el_angle_offset.append(foc.el_angle_offset)
        self.torque_limit.append(foc.torque_limit)

        return (cmd_a, cmd_b, cmd_c)
