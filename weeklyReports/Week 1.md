# Week 1 Report (10th June-14th June)

## Literature review

### Different types of Fire suppression systems under consideration

#### 1. SAPPHIRE® PLUS 70-BAR SYSTEM

[  https://www.ansul.com/en/us/pages/ProductSeries.aspx?ProductType=Clean+Agent+Systems ]

Chemical used: 3M Novec 1230 [Refer 3M pdf for detailed properties] or dodecafluoro-2-methylpentan-3-one
Pressure: 70bar (mostly 25-42 bar pressurized with help of nitrogen) [ref: https://www.tycoifs.co.uk/how-we-can-help/protect-your-business/fire-suppression-systems-marine/gaseous-systems-sapphire-novec-1230/ ]
Container size range: 8ltr to 180 ltr
Release time: within 10sec (depends on combination of pressure and valve)
Method: Fully flooding, mainly absorbs heat and has some chemical interference with the flame. Upon discharge forms a gaseous mixture with air
Specialties:
1.	Higher fill densities reduces no. of cylinders
2.	Increased pressure helps in keeping containers far away and selector valves possibles
3.	pipe size may also be reduced to save money and design flexibility
4.	good for class A,B and C fire hazard
5.	Evaporates 50 times faster than water
6.	Electrically non-conductive
7.	Very Low water solubility
8.	principal atmospheric sink for Novec 1230 agent is photolysis

Concentrations: 4.5 – 5.9% (safe up to 10%)
ODP=0, ALT= 3-5 DAYS, GWP=1.0

Safety: 
NOVEC 1230 decomposes at temperatures in excess of 500°C and it is therefore important to avoid applications involving hazards where continuously hot surfaces are involved. Upon exposure to the flame, NOVEC 1230 will decompose to form halogen acids. Their presence will be readily detected by a sharp, pungent odor before maximum exposure levels are reached.

Each container assembly comprises of a container, valve, siphon tube, safety fittings and container nameplate. Containers are painted signal red to identify them as fire protection equipment.
The valve is designed for high discharge rates to allow the full container contents to be released within ten seconds.

Thermal properties:
![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/Thermal%20properties%20of%20Novec.png "Thermal Properties")


Variation in vapor pressure and density with temperature:

![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/Novec%20Vapor%20pressure%20vs%20temp.png "Vapor pressure vs temp")
![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/Novec%20Density%20vs%20Temp.png "Density vs Temp")

Even though liquid has enough vapor concentration to suppress fire. At 39% volume concentration in air it reaches saturation. In applications concentration used is less than 6%


Novec 1230 fluid is such that it would support an extinguishing concentration of 5vol% at a temperature as low as –16°C (3°F). Water does not support a 5vol% concentration in air until the temperature exceeds 33°C(91°F). This can be seen in following graph.

![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/Novec's%20Saturation%20comparison%20with%20water.png "Saturation comparison")



https://instor.com/blog/which-data-center-fire-suppression-solution-works-best-for-you/

#### Inergen Systems
Inergen is a breathable gas comprised of 52% Nitrogen, 40%Argon and 8% Carbon Dioxide. It suppresses a fire by lowering the oxygen level in a room below the level needed to support combustion, but still leaves the oxygen level high enough to sustain life, while suppressing the fire. The release time for Inergen is 60 seconds, which is the longest release time of these three options. Inergen systems require more intricate pipe systems that also require more pressure to be maintained within the system.

#### FM-200 Systems
FM-200 systems also suppress fires by removing heat from the room. They can be used with people in the room with few adverse effects. FM-200 systems are also much quicker than Inergen systems to reach full fire suppression levels.  FM-200 has a much better potential for having a global warming effect than the other two systems.
All three of these different clean agent systems can be used while employees remain inside of the room, but Inergen seems to have more of a risk factor for some people because it lowers the oxygen level, which may not be ideal for all people. All three of these systems require floor space and are costlier to install and maintain than water systems. In addition to these disadvantages, clean agent systems can be linked to the EPO (Emergency Power Off) System. This was a requirement of the NFPA until 2011.

#### CO2 Agent: 
Carbon dioxide is an effective fire suppressing agent that can be used on many types of fires. It is effective for surface fires, such as flammable liquids and most solid combustible materials. It expands at a ratio of 450 to 1 by volume. For fire suppression purposes, the discharge is designed to raise the carbon dioxide concentration in the hazard. This displaces the air, which contains oxygen that supports combustion, and results in fire suppression. Other attributes are its high degree of effectiveness, its excellent thermal stability, and its freedom from deterioration. It is electrically non-conductive, and leaves no residue to clean up after discharge.

### DIFFERENT MODEL USED FOR ANALYSIS

#### Modeling of the Flow Properties and Discharge of Halon Replacement Agents (by P.J. DiNenno, E.W. Forsseli, M.J. Ferreira, C.P. Hanauska, and B.A. Johnson, Hughes Associates)

##### * Assumptions

1.Conditions in cylinder (P, T, composition) = f( initial conditions, outage fraction)

  where outage fraction= fraction of initial mass left; This assumption ignores impact of K.E. of fluid leaving cylinder on     the cylinder energy balance
  
2.Quasi-steady Flow

  (average flow rate over a small time step is equal to the flow rate that would exist if the cylinder conditions were held steady during the time step)
  
3.Flow through pipe is homogeneouos (both phases at same velocity with one phase totally dispersed in other)

4.Heat transfer through the pipes is considered insignificant. (This can't be done for our model)
  
##### * Method used (HFLOW with small extension for application)

1.Flow rate estimated approximately (guess)

2.Network is stepped through to determine condition at nozzle

3.Use Energy and momentum balance eqn to find P and T at each node.
(branches in the network stepped by 34.5KPa pressure drops and distance travelled determined through momentum balance but some adjustment can be done in the pressure to keep the distance between nodes less than 8cm. The distance aims at easy identification of location for sonic condition  )

4.Estimated flow rate is then refined by comparison to deteremined flow through the nozzles. Iterations done.

5.Elapsed time measured using mass balance (Mass that has already left cylincer/ mass flow rate)

##### *Variation in pressure vs Time

![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/Discharge%20flow%20phenomenon.png "Discharge flow phenomenon")


In Normal Scenario (a):

Start – (1): agent not reached the nozzle yet, progressing through network at sonic conditions

(1)-(2): After the pipe gets full at (1), pressure and mass in the pipe network starts building up to reach peak value. (In this section flow rate through the network starts dropping whereas flow from the nozzle starts increasing)

(2)-(3): Mass flow rate out of cylinder is equal to that of nozzle.

(3)-(4): Vapor front from the cylinder starts flowing in the pipe network in the same way as the liquid front started initially. The low flow rate of liquid ahead of the vapor in the network controls and prevents vapor from increasing its flow rate.  

After(4): Once nozzle is cleared of liquid, the flow rate of vapor from cylinder increases rapidly and then falls off (flow rate from cylinder = nozzle)

Other scenarios:

(b).Cylinder runs out of liquid prior to reaching peak pressure in pipe:
In such case the pressure at nozzle will continue to rise but the apex won’t be that high. Also steady flow section never comes i.e (2)-(3)

(c).Nozzle does not control the flow rate or controls it completely: 
In extreme cases pipe peak pressure already occurs before liquid front reaches nozzles. i.e. section (1)-(2) is skipped.

##### *Flow chart/Algorithm of the model

![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/Flow%20chart%20algorithm.PNG "Flow chart algorithm")



##### *Accuracy and Limitations of Model

1.Model can measure the discharge time based on nozzle liquid runout with reasonable accuracy and  its easy to modify for halon kind of gases.

2.Model halts execution if the cylinder runs out of liquid prior to liquid reaching the last nozzle. This limits the pipe volume roughly  to an NFPA 12A 70% agent in pipe. (pipe volume/ initial cylincder liquid volume = 1.6)

3.Needs to tested on more orientation of pipe networks. (here it is tested against 2 Tee orientation: Horizontal bull-head and horizontal side flow)

4.Largest flow splits that the model can handle is also not known for different pipe networks. (The two orienations used here have maximum flow split around 90%/10%)


##### *Original HFLOW Model (refered in this document) is taken from the following research paper: Flow of Nitrogen- Pressurized Halon 1301 in Fire Extinguishing syatems by Elliot et. al.


![alt text](https://github.com/Yashwantyogi/twoPhaseSimulationInternProject/blob/master/Images%20and%20Graph/HFLOW%20model.PNG "Original HFLOW model")
Note: The numbers in the diagram show the sequence of steps the program goes through

##### Will add more details of this research paper and its model in next Weekly report
