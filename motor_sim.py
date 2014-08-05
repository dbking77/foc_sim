#!/usr/bin/env python
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


"""
Simulates sinusoidal 3-phase brushless motor with using very basic priciples
and treats phase inductances separately and does not take any Q-D shortcuts
doesn't currently support sailent motors or non-sinusoidal motors.

The purpose of this simulator is to allow testing and debugging of different
control algoritms and control gains on different motor/motor parameters.
The simulator combinded with ploting tools should allow easier insight into
motor and controller operation that can be had from running on hardware system

Simulator uses first-order (Euler or RK1) method to simulate motor differential
equations.  The Euler method needs a smaller simulation timestep (than other methods)
is need to make simulation stable.  However, the Euler method is much easier to implement.

There are two time-constants to worry about for simulator.
The motor inductance time-constant, and the motor inertial time-constant.
In practice most real motors have an indutance time-constant that is always smaller
than motor inertial time-constant so it can be ignored in this case.

The current timestep is good for most types of motors however might need to be reduced
if motors has very small inductance (< 100uH).

While simulator is meant to be somewhat generic, it does not implement
all possible modes of operation as configuration parameters, instead modify the
simulator code to support type of motor or controller operation that is of interest.
"""


from math import sin, cos, pi, sqrt
import pylab
from foc import FOC_Controller, FOC_CWrapper, alignmentRatio

################################################################################
# Block-Commutation Versus simulation
# Parameters of maxon-EC-Motors Driven By Servo Amplifiers With Sinus-Commutation
#
# The torque equation we use for simulation
#   T = Kt * (Ia*sin(A) + Ib*sin(A-2*pi/3) + Ic*sin(A+2*pi/3))
#
# In this equation Kt is the torque constant.
# Unfortunately this is not the torque constant that Maxon provides (Kb).
#
# According to Maxon, the torque constant Kb is determine using the average torque
# average torque produced when using block communtation.
#
# In block commutation, two phases are conducting at a time, while the third phase
# has zero current.  Each motor electrical revolution consists of 6-steps. 
# The phases are energized so that they maximize torque production during that step
# of the motor revolution.
#
# During a single commutation step (1/6th of a rotation) the currents would be constant
# For one of the block commutation steps Ia = -Ib and Ic = 0
# If commutation is done correctly, these currents should occur when angle (A) is between
# 30 and 90 degrees (pi/6 and pi/2)
#
# The average torque produced over the 1/6th (1/6 revolution == pi/3 radians) will be
#   Tavg = 1/(pi/3) * Kt * integral(Ia*sin(A) + Ib*sin(A-2*pi/3) + Ic*sin(A+2*pi/3), A, pi/6, pi/2)
# 
# For block commutation Ia=-Ib, Ic = 0
#   Tavg = 1/(pi/3) * Kt * integral(Ia*sin(A) - Ia*sin(A-2*pi/3), dA, pi/6, pi/2)
#
# Intrgral is easy to solve
#   Tavg = 1/(pi/3) * Kt * (-Ia*cos(A) + Ia*cos(A-2*pi/3), | A(pi/6, pi/2)
#   Tavg = 1/(pi/3) * Kt * (-Ia*cos(pi/2) + Ia*cos(pi/2-2*pi/3)) - (-Ia*cos(pi/6) + Ia*cos(pi/6-2*pi/3))
#   Tavg = 1/(pi/3) * Kt * (-Ia*0 + Ia*cos(pi/6)) - (-Ia*cos(pi/6) + Ia*0)
#   Tavg = 1/(pi/3) * Kt * 2*Ia*cos(pi/6)
#
# cos(pi/6) = sqrt(3)/2
#   Tavg = 1/(pi/3) * Kt * 2*Ia*sqrt(3)/2
#   Tavg = 1/(pi/3) * Kt * Ia*sqrt(3)
#   Tavg = 3/pi * Kt * Ia * sqrt(3)
#   Tavg = 3/pi * Kt * sqrt(3) * Ia
#
#
# If Tavg = Kb * Ia then:
#   Kb * Ia = Tavg = 3/pi * Kt * sqrt(3) * Ia
#   Kt = Kb * pi/3 / sqrt(3)
#
# Conclusion : The Kt value used in the simulated torque equation should be
#   Kt = Kb * pi/3 / sqrt(3)
#
################################################################################


################################################################################
# Why subtract out average voltage?
#   va = dIa/dt*L + Ia*R + backemf_a + v_center
#   vb = dIb/dt*L + Ib*R + backemf_b + v_center
#   vc = dIc/dt*L + Ic*R + backemf_c + v_center
#   Ia + Ib + Ic = 0
#
#  Also
#   dIa + dIb + dIc = 0
#   dIa = -dIb - dIc
#
#
#  v_center = va - (-dIb-dIc)/dt*L - Ia*R - backemf_a
#  v_center = va + (dIb+dIc)/dt*L - (-Ib-Ic)*R - backemf_a
#  v_center = va + dIb/dt*L + dIc/dt*L + Ib*R + Ic*R - backemf_a
#  v_center = va + (dIb/dt*L + Ib*R) + (dIc/dt*L + Ic*R) - backemf_a
#
#  vb - backemf_b - v_center = dIb/dt*L + Ib*R
#  vc - backemf_c - v_center = dIc/dt*L + Ic*R  
#
#  v_center = va + vb - backemf_b - v_center + vc - backemf_c - v_centerx - backemf_a
#  3 * v_center = (va-backemf_a) + (vb-backemf_b) + (vc-backemf_c
#  v_center = [(va-backemf_a) + (vb-backemf_b) + (vc-backemf_c] / 3
#  v_center = [(va+vb+vc) - (backemf_a + vb-backemf_b + vc-backemf_c)] / 3
#
#  since backemf are 3-phased they sum to 0, so finally
#    v_center = (va + vb + vc) / 3
#
################################################################################





# Parameters for Maxon EC-Flat 45 30-Watt 24V
# Parameters are taken from Maxon data-sheet, however these parameters
# are based on block commutation control of motors
class Maxon_339283:
    def __init__(self):
	self.torque_constant = 51e-3 # Nm/A == V/(rad/s)
	self.terminal_resistance = 5.03
	self.terminal_inductance = 2.24e-3 #mH
	self.rotor_inertia = 92.5 * 1e-7  # gcm2 = kg/1000 * (m/100 * m/100) = kg/m2/1e7

	#no load speed = 4380 rpm
	#no load speed = 458.67 rad/s
	#no load current = 73mA
	#no load torque = 73mA * 51  = 3.723e-3 Nm
	self.friction = 3.723e-3 / 458.67  # Nm/(rad/s)

        self.max_torque = 54.7e-3

        self.pole_pairs = 8



# These controller gains work well for 17.5kHz control loop rate and Maxon_339283
class Maxon_339283_Controller_Gains:
    def __init__(self):
        # d_gains for direct current controller,
        # q_gains for quadrature current controller
        # first value is proportial term for PI controller
        # second value is integral term for PI controller
        self.d_gains = (0.4, 800.)
        self.q_gains = (0.4, 800.)





# Simulates 3-phase brushless motor with sinusoidal back-EMF using parameters from motor class.
# Assume motor phases are Wye (or star) connected so are no circulating currents.
# Currently, cannot simulate field-weaking, relectance torque, or no-sinusoid back-EMF
# Keeps arrays of motor states, inputs and outputs to allow easy plotting of motor operation later,
# However these arrays can easily get too big if simulation is run too long.
class MotorSim:
    def __init__(self, motor, load_inertia, start_position=0.0):
	self.motor =motor #motor parameter (backemf contant, torque constant, terminal resistance, terminal inductance, inertia, etc..)
        self.load_inertia = load_inertia
	self.ia = [0.0]  # phase currents
	self.ib = [0.0]
	self.ic = [0.0]
	self.vxa = [0.0]  # terminal voltages
	self.vxb = [0.0]
	self.vxc = [0.0]
	self.va = [0.0]  # phase voltages
	self.vb = [0.0]
	self.vc = [0.0]
        self.backemf_a = [0.0]  # back-emf voltages
        self.backemf_b = [0.0]
        self.backemf_c = [0.0]
	self.torque = [0.0]   # torque produced  
	self.velocity = [0.0] # motor velocity (radians/sec)
	self.position = [start_position] # motor position (radians)
	self.t = [0.0]
    
    def update(self, t, dt, vxa, vxb, vxc):
	""" Updates motor state (phase currents, velocity, position) given excited voltages on each phase
	t is the current time
	dt is the timestep
	"""
	motor = self.motor
	ia = self.ia[-1]
	ib = self.ib[-1]
	ic = self.ic[-1]
	
	# phase voltages are not absolute (since there is not ground)
	# as such sum of voltages should be 0.
	# to acomplish this, subtrace out average
	vcenter = ( vxa + vxb + vxc ) / 3.0
	va = vxa - vcenter
	vb = vxb - vcenter
	vc = vxc - vcenter

	# resistance of a single phase winding
	R = motor.terminal_resistance * 0.5	

	# determine voltages that are seen by inductances of each phase
	# by removing resitance and back-emf voltages
	Lva = va - self.backemf_a[-1] - ia * R
	Lvb = vb - self.backemf_b[-1] - ib * R
	Lvc = vc - self.backemf_c[-1] - ic * R

	# inductance of a single phase
	L = motor.terminal_inductance * 0.5

	# volatage will induce current change in phase current
	ia += Lva / L * dt
	ib += Lvb / L * dt
	ic += Lvc / L * dt

	# sum of current should be 0, kill non-zero (circulating) currents
	# this assumes the motor is Wye-wound and not Delta wound
        iavg = (ia+ib+ic)/3.0
        ia -= iavg
        ib -= iavg
        ic -= iavg

	# sin of motor angle will be important for a couple calcs
        el_angle = self.position[-1] * self.motor.pole_pairs
	sa = sin(el_angle)
	sb = sin(el_angle-pi*2/3)
	sc = sin(el_angle+pi*2/3)

	# determine torque
	# note that the torque_constant is based on block commutation.
	# read notes about motor torque with sinusoidal control to
	# get a better idea scale constant comes from
	scale  = 0.60459978807807258  # pi/3 / sqrt(3)
	torque_constant = motor.torque_constant * scale
	torque = (ia*sa + ib*sb + ic*sc) * torque_constant

	# determine new motor velocity (include viscous friction)
        velocity = self.velocity[-1]
        velocity += (torque - velocity*motor.friction)/(motor.rotor_inertia+self.load_inertia)*dt

	# determine new motor position
	position = self.position[-1] + velocity*dt

	# determine new back-emf voltages (for next time)
	# for motor V/(rad/s) == A/Nm == 1/(Nm/A)
	#
	# Nm/A = Nm/A * s/s
	# Nm/A = ((Nm/s)*s)/A
	# Nm/A = (W*s)/A
	# Nm/A = (W/A)*s
	# Nm/A = (V)*s
	# Nm/A = V/(1/s)
	# Nm/A = V/(rad/s)
	voltage_constant = torque_constant
	backemf_a = sa * velocity * voltage_constant
	backemf_b = sb * velocity * voltage_constant
	backemf_c = sc * velocity * voltage_constant

	# record 
	self.vxa.append(vxa)
	self.vxb.append(vxb)
	self.vxc.append(vxc)
	self.va.append(va)
	self.vb.append(vb)
	self.vc.append(vc)
	self.ia.append(ia)
	self.ib.append(ib)
	self.ic.append(ic)
	self.backemf_a.append(backemf_a)
	self.backemf_b.append(backemf_b)
	self.backemf_c.append(backemf_c)
	self.torque.append(torque)
	self.velocity.append(velocity)
	self.position.append(position)
	self.t.append(t)

        return (position, velocity, ia, ib, ic)



def sim():
    # Motor parameters can be switched here pretty quickly
    motor = Maxon_339283()
    load_inertia = 0.0 * 1e-7 #0.0 # Extra intertial load to use when simulating motor (kg*m^2)

    gains = Maxon_339283_Controller_Gains()
    d_gains = gains.d_gains
    q_gains = gains.q_gains

    alignmentRatio(motor)

    # motor simulator 
    motor_sim = MotorSim(motor, load_inertia)

    
    # Simulate starting offset to absolte motor position measurement 
    #  (for testing incremental alignment)
    el_angle_offset = 0.0 #-20.0 * pi/180.
    angle_offset = el_angle_offset / motor.pole_pairs

    # Supply voltage, controller outputs PWM duty values that should be in range 0.0-1.0.
    # the simulator multiple the supply voltage to simulate a 3-phase inverter
    supply_voltage = 24.0

    # For FOC controller this is the initial target torque that is passed in.
    # if more complex torque profile or higher level position / velocity controller is
    # needed, implement it in simulatin loop
    target_torque = 0.65 * motor.max_torque

    # Sim time
    t = 0.0

    # Inverter PWM frequency. Simulator does not simulate actual PWM (on-off with specific duty cycle)
    # However simulator does change output voltage in discrete times steps to mirror actual PWM based
    # inverter that can only change its output at discrete time-steps.
    pwm_dt = 1.0 / 35e3 

    # How many PWM cycles occur before FOC update is run.
    # If this value is 1 the FOC is recalculated after every PWM cycle
    # This code makes assumption that FOC is run synchronously to PWM cycle.
    #
    # However, in some controllers the FOC is run asynchornously to PWM updates,
    # This decouples PWM frequency from FOC update frequency, but at cost
    # of less repeatable performance since FOC and PWM update are sliding in and
    # out of phase.
    #
    # To simulate a FOC update that is asynchrounous to PWM, change code in
    # main simulation loop to keep separate update time for FOC controller.
    foc_update_cycles = 1
    foc_dt = pwm_dt * foc_update_cycles

    # Some FOC controllers enforce supply currnet (or supply power limit)
    supply_current_limit = 1.2 # negative value = no current limit

    foc = FOC_Controller(foc_dt, motor, d_gains, q_gains, supply_current_limit)
    #foc = FOC_CWrapper(foc_dt, motor, d_gains, q_gains, supply_current_limit)
    foc.reset(motor_sim.velocity[-1], supply_voltage) 

    #
    pwm_t = []
    pwm_cmd_a = []
    pwm_cmd_b = []
    pwm_cmd_c = []

    if True:
        cmd_a,cmd_b,cmd_c = foc.update(t, 0.0, motor_sim.position[-1], motor_sim.velocity[-1], 0., 0., 0., supply_voltage)
        vxa = supply_voltage*cmd_a
        vxb = supply_voltage*cmd_b
        vxc = supply_voltage*cmd_c
    else:
        vxa,vxb,vxc = (0.0, 0.0, 0.0)
        cmd_a,cmd_b,cmd_c = (0.0, 0.0, 0.0)

    motor_sim_substeps = 8
    motor_sim_dt = pwm_dt/motor_sim_substeps
    runtime = 15e-3 
    while t < runtime:
        for j in range(foc_update_cycles):
            for i in range(motor_sim_substeps):
                position,velocity,ia,ib,ic = motor_sim.update(t, motor_sim_dt, vxa, vxb, vxc)
                t+=motor_sim_dt
            # calulate command voltages here to simulate 1 cycle control loop delay
            vxa = supply_voltage*cmd_a
            vxb = supply_voltage*cmd_b
            vxc = supply_voltage*cmd_c
            pwm_t.append(t)
            pwm_cmd_a.append(cmd_a)
            pwm_cmd_b.append(cmd_b)
            pwm_cmd_c.append(cmd_c)

        cmd_a,cmd_b,cmd_c = foc.update(t, target_torque, position+angle_offset, velocity, ia,ib,ic, supply_voltage)


    ts = pylab.array(motor_sim.t) * 1e3
    tc = pylab.array(foc.t) * 1e3

    va = pylab.array(motor_sim.va)
    vb = pylab.array(motor_sim.vb)    
    vc = pylab.array(motor_sim.vc)
    backemf_a = pylab.array(motor_sim.backemf_a)
    backemf_b = pylab.array(motor_sim.backemf_b)
    backemf_c = pylab.array(motor_sim.backemf_c)
    ia = pylab.array(motor_sim.ia)
    ib = pylab.array(motor_sim.ib)
    ic = pylab.array(motor_sim.ic)
    vxa = pylab.array(motor_sim.vxa)
    vxb = pylab.array(motor_sim.vxb)    
    vxc = pylab.array(motor_sim.vxc)

    cmd_d = pylab.array(foc.cmd_d)
    cmd_q = pylab.array(foc.cmd_q)
    id = pylab.array(foc.measured_id)
    iq = pylab.array(foc.measured_iq)

    pwm_t = pylab.array(pwm_t) * 1e3
    pwm_cmd_a = pylab.array(pwm_cmd_a)
    pwm_cmd_b = pylab.array(pwm_cmd_b)
    pwm_cmd_c = pylab.array(pwm_cmd_c)

    if True: #FOC stuff
        pylab.figure("FOC")
        pylab.subplot(3,1,1)
        pylab.plot(tc, foc.cmd_a, 'r.-', label='cmd a')
        pylab.plot(tc, foc.cmd_b, 'g.-', label='cmd b')
        pylab.plot(tc, foc.cmd_c, 'b.-', label='cmd c')

        pylab.plot(pwm_t, pwm_cmd_a, 'r.--', label='pwm cmd a')
        pylab.plot(pwm_t, pwm_cmd_b, 'g.--', label='pwm_cmd b')
        pylab.plot(pwm_t, pwm_cmd_c, 'b.--', label='pwm_cmd c')

        pylab.legend()
	pylab.title(foc.name)
        pylab.xlabel('time ms')    
        pylab.subplot(3,1,2)
        pylab.plot(tc, 1e3*pylab.array(foc.target_id), 'r', label='target Id (mA)')
        pylab.plot(tc, 1e3*pylab.array(foc.measured_id), 'g', label='measured Id (mA)')
        pylab.legend()
        pylab.xlabel('time ms')    
        pylab.ylabel('mAmps')
        pylab.subplot(3,1,3)
        pylab.plot(tc, foc.target_iq, 'r', label='target Iq')
        pylab.plot(tc, foc.measured_iq, 'g', label='measured Iq')
        pylab.legend()
        pylab.xlabel('time ms')    
        pylab.ylabel('Amps')    

    
    if False : # PI controllers
        pylab.figure("PI Controllers")
        d_ctrl = foc.d_ctrl
        d_ctrl_t = pylab.array(d_ctrl.t) * 1e3
        pylab.subplot(2,1,1)
        pylab.plot(d_ctrl_t, d_ctrl.i_term, 'r', label='d_ctrl : i-term')
        pylab.plot(d_ctrl_t, d_ctrl.output, 'b', label='d_ctrl : output')
        pylab.legend()
        pylab.xlabel('time ms')    
        q_ctrl = foc.q_ctrl
        q_ctrl_t = pylab.array(q_ctrl.t) * 1e3
        pylab.subplot(2,1,2)
        pylab.plot(q_ctrl_t, q_ctrl.i_term, 'r', label='q_ctrl : i-term')
        pylab.plot(q_ctrl_t, q_ctrl.output, 'b', label='q_ctrl : output')
        pylab.plot(tc, foc.cmd_q_filt, 'g', label='cmd_q_filt')

        pylab.legend()
        pylab.xlabel('time ms')    

    print('peak current', max(ia.max(), ib.max(), ic.max()))

    pylab.figure("Sim Currents, Voltage, and Power")
    pylab.subplot(3,1,1)
    pylab.plot(ts, ia, 'r', label='Ia')
    pylab.plot(ts, ib, 'g', label='Ib')
    pylab.plot(ts, ic, 'b', label='Ic')
    pylab.legend()
    pylab.xlabel('time ms')

    pylab.subplot(3,1,2)
    pylab.plot(ts, backemf_a, 'r', label='backemf A')
    pylab.plot(ts, backemf_b, 'g', label='backemf B')
    pylab.plot(ts, backemf_c, 'b', label='backemf C')
    pylab.legend()


    R = motor.terminal_resistance*0.5
    rpower = R*(ia*ia+ib*ib+ic*ic)

    # find inductive power by taking delta on inductor stored energy
    # E = 1/2*L*(I**2)
    L = motor.terminal_inductance*0.5
    Lenergy = 0.5*L*(ia*ia+ib*ib+ic*ic)
    Lpower = pylab.diff(Lenergy)
    Lpower = pylab.append([0.0],Lpower)/motor_sim_dt


    powerA = va*ia
    powerB = vb*ib
    powerC = vc*ic
    elec_power = powerA + powerB + powerC

    powerQD = (id*cmd_d + iq*cmd_q)*supply_voltage / sqrt(2.0)

    torque   = pylab.array(motor_sim.torque)
    velocity = pylab.array(motor_sim.velocity)
    mech_power = torque * velocity

    
    pylab.subplot(3,1,3)
    pylab.plot(ts, mech_power, 'r', label='mechanical power')
    pylab.plot(ts, elec_power, 'g', label='electical power')
    pylab.plot(tc, powerQD, '0.5', label='QD power')
    pylab.plot(ts, rpower, 'b', label='heating power')
    pylab.plot(ts, Lpower, 'c', label='inductive power')
    pylab.plot(ts, rpower+mech_power+Lpower, 'k--', label='heating + mech + inductive')
    #pylab.plot(ts, powerA, 'r', label='power A')	
    #pylab.plot(ts, powerB, 'g', label='power B')	
    #pylab.plot(ts, powerC, 'b', label='power C')	    
    pylab.legend()


    print("power ratio ", (powerQD[-1]/elec_power[-1]))

    pylab.figure("Torque, Velocity, and Position")
    pylab.subplot(3,1,1)
    pylab.plot(tc, pylab.array(foc.input_target_torque) * 1e3, 'c', label='input target torque')
    pylab.plot(tc, pylab.array(foc.target_torque) * 1e3, 'g', label='FOC target torque mNm')
    pylab.plot(ts, torque * 1e3, 'r--', label='motor sim torque mNm')
    pylab.plot(tc, pylab.array(foc.measured_torque) * 1e3, 'b', label='FOC measured torque mNm')
    ylim = pylab.ylim()
    pylab.plot(tc, pylab.array(foc.torque_limit) * 1e3, 'k--', label='dynamic torque limit mNm')
    # torque limit can be really large in some cases, don't allow plot to zoom out because of this
    yrange = ylim[1]-ylim[0]
    pylab.ylim( (ylim[0]-yrange*0.1, ylim[1]+yrange*0.1)  )


    pylab.ylabel("mNm")
    pylab.legend()
    pylab.subplot(3,1,2)
    pylab.plot(ts, velocity*60/(2*pi), 'r', label='velocity (RPM)')
    pylab.legend()
    pylab.subplot(3,1,3)
    pylab.plot(ts, pylab.array(motor_sim.position)*180./pi, 'r', label='position (degrees)')
    pylab.legend()


    # EL Angle offset error estimate
    if True:
        pylab.figure("Angle Offset")
        pylab.subplot(3,1,1)
        pylab.plot(tc, foc.error_vd, 'r', label='Vd error')
        pylab.plot(tc, foc.cmd_vd, 'g', label='Vd cmd')
        pylab.plot(tc, foc.est_vd, 'b', label='Vd est')
        pylab.legend()
        pylab.subplot(3,1,2)
        pylab.plot(tc, foc.backemf_voltage, 'r', label='backemf voltage')
        pylab.plot(tc, foc.est_vd, 'b', label='Vd est')
        pylab.legend()
        pylab.subplot(3,1,3)
        pylab.plot(tc, pylab.array(foc.est_el_angle_offset)*180/pi, 'r', label='est_angle_offset (degrees)')
        pylab.plot(tc, pylab.array(foc.el_angle_offset)*180/pi, 'g', label='angle_offset (degrees)')
	pylab.plot(tc, (el_angle_offset-pylab.array(foc.el_angle_offset))*180/pi, 'k--', label='actual - angle_offset (degrees)')
        pylab.plot([0,tc[-1]], [el_angle_offset*180/pi, el_angle_offset*180/pi], 'b-*', label='actual el angle offset (degrees)')
        pylab.legend()
        print foc.measured_iq[-1]


    if False:
        pylab.figure()
        current_vs_torque = ( ((ia*ia + ib*ib + ic*ic)*0.5)**0.5 ) / torque * motor.torque_constant
        pylab.plot(ts, current_vs_torque, label="Current/Torque sqrt(Ia^2+Ib^2+Ic^2)/(2*torque)")
        pylab.legend()
        print "Current/Torque sqrt(Ia^2+Ib^2+Ic^2)/(2*torque) : ", current_vs_torque[20:].mean()

    pylab.show()



if __name__ == "__main__":
    sim()

