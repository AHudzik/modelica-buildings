within Buildings.HeatTransfer.Conduction;
model SingleLayerHM

  Modelica.SIunits.MassFlowRate m_flow[nSta + 1];
  Modelica.SIunits.HeatFlowRate Q_flow[nSta + 1];
  Modelica.SIunits.Temperature T[nSta](each start=T_ini)
    "Temperature at the states";

  Buildings.HeatTransfer.Types.VapourDiffusionCoefficient kd[
    nSta + 1];

  Modelica.SIunits.Pressure pw[nSta] "Vapour Pressure at the state";

  Buildings.HeatTransfer.Types.RelativeHumidity phi[nSta](
    each start=phi_ini,
    max=1.0,
    min=0.0) "Relative humidity";

  Buildings.HeatTransfer.Types.MoistureStorageCapacity dw_dphi[
    nSta] "Derivative of water content over relative humidity.";

  Modelica.SIunits.Pressure psat[nSta]
    "equilibrum vapor pressure at the states";

  Modelica.SIunits.MassFraction Xi_outflow[nSta] "mass fraction at the state";

  Modelica.SIunits.MassConcentration w[nSta] "Water content at the states";
                                            //(each start=w_ini)

  Modelica.SIunits.HeatCapacity Cm[nSta]
    "Heat capacity of moist material at the states";

  Modelica.SIunits.ThermalConductivity lambda[nSta];

  Modelica.SIunits.ThermalConductance UA[nSta + 1];

  Buildings.HeatTransfer.Types.LiquidConductionCoefficient
    D_phi[nSta + 1] "liquid conduction coefficient [kg/ms]";

  parameter Modelica.SIunits.Area A=1 "Heat transfer area";

  parameter Boolean activatesuction;

  parameter Buildings.HeatTransfer.Types.RelativeHumidity
    phi_ini=0.0 "Initial value of relative humidity at the nodes";
  parameter Modelica.SIunits.Temperature T_ini=293.15;

  parameter Integer switch_w "switch for the water content"
    annotation (Dialog(tab="HM"));
  parameter Integer switch_lamb "switch for the heat conductivity"
    annotation (Dialog(tab="HM"));
  parameter Integer switch_dw
    "switch for the liquid transport coefficient for suction "
    annotation (Dialog(tab="HM"));

  replaceable parameter Data.BaseClasses.HygroThermalMaterial material
    "Material from Data.Solids, Data.SolidsPCM or Data.Resistances" annotation (
    Evaluate=true,
    choicesAllMatching=true,
    Placement(transformation(extent={{60,60},{80,80}})));

protected
  Modelica.SIunits.TemperatureDifference dT "port_a.T - port_b.T";

  final parameter Modelica.SIunits.Volume V=A*material.x/nSta;
  final parameter Integer nSta(min=1) = material.nSta
    "Number of state variables";

  parameter Modelica.SIunits.Mass m=A*material.x*material.d/material.nSta
    "Mass associated with the temperature state";
  parameter Modelica.SIunits.HeatCapacity Cd=m*material.c
    "Heat capacity of dry material associated with the temperature state";
  parameter Real b=0.8*(material.w_80 - material.w_f)/(material.w_80 - 0.8*
      material.w_f);
  parameter Modelica.SIunits.Length d_node = material.x / nSta
    "thickness of numerical node";

  constant Modelica.SIunits.Pressure patm=101325 "Atmospheric pressure";
  constant Modelica.SIunits.SpecificHeatCapacity Cw=4185.0;
  constant Modelica.SIunits.Enthalpy hv=2500000
    "Evaporation enthalpy of water [J/Kg]";

public
  Buildings.HeatTransfer.Interfaces.HeatMassPort_b heatMassPort_b
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  Interfaces.HeatMassPort_a heatMassPort_a
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

equation
  dT = heatMassPort_a.heatPort.T - heatMassPort_b.heatPort.T;

  heatMassPort_a.heatPort.Q_flow = +Q_flow[1];
  heatMassPort_a.heatPort.T - T[1] = Q_flow[1]/UA[1];
  heatMassPort_b.heatPort.Q_flow = -Q_flow[nSta + 1];
  T[nSta] - heatMassPort_b.heatPort.T = Q_flow[nSta + 1]/UA[nSta + 1];

  heatMassPort_a.massPort.m_flow = +m_flow[1];
  Buildings.Utilities.Psychrometrics.Functions.pW_X(heatMassPort_a.massPort.Xi_outflow)
     - pw[1] = m_flow[1]/(kd[1]*A);
  heatMassPort_b.massPort.m_flow = -m_flow[nSta + 1];
  pw[nSta] - Buildings.Utilities.Psychrometrics.Functions.pW_X(heatMassPort_b.massPort.Xi_outflow)
    = m_flow[nSta + 1]/(kd[nSta + 1]*A);

  for i in 1:nSta loop
    dw_dphi[i] = Buildings.HeatTransfer.Conduction.Functions.dw_dphi(
      phi[i],
      material.w_f,
      b,
      material.sorp_tab,
      material.por,
      switch_w);

    lambda[i] = Buildings.HeatTransfer.Conduction.Functions.lambdaHM(
      material.d,
      material.k,
      w[i],
      material.b_h,
      material.lamb_tab,
      switch_lamb);

    w[i] = Buildings.HeatTransfer.Conduction.Functions.sorption(
      phi[i],
      material.w_f,
      b,
      material.sorp_tab,
      material.por,
      switch_w);

    Xi_outflow[i] = Buildings.Utilities.Psychrometrics.Functions.pW_X(pw[i]);

  end for;

  // WATER VAPOUR DIFFUSION COEFFICIENT
  kd[1] = Buildings.HeatTransfer.Conduction.Functions.delta_L(T[1])/(d_node/2)/material.mu;
  for i in 2:nSta loop
    kd[i] = 1.0/ ( 1.0/(
      Buildings.HeatTransfer.Conduction.Functions.delta_L(T[i])/(d_node/2)/material.mu
       +
      Buildings.HeatTransfer.Conduction.Functions.delta_L(T[i - 1])/(d_node/2)/material.mu));
  end for;
  kd[nSta + 1] = Buildings.HeatTransfer.Conduction.Functions.delta_L(T[nSta])/(d_node/2)/material.mu;

  // LIQUID WATER DIFFUSION COEFFIIENT
  D_phi[1] = dw_dphi[1]*Buildings.HeatTransfer.Conduction.Functions.Dw(
    w[1],
    material.w_f,
    material.A_layer,
    material.dww_tab,
    material.dws_tab,
    switch_dw,
    activatesuction)/(d_node/2);
  for i in 2:nSta loop
    D_phi[i] = 1.0/((d_node/2)/(dw_dphi[i]*
      Buildings.HeatTransfer.Conduction.Functions.Dw(
      w[i],
      material.w_f,
      material.A_layer,
      material.dww_tab,
      material.dws_tab,
      switch_dw,
      activatesuction)) + (d_node/2)/(dw_dphi[i - 1]*
      Buildings.HeatTransfer.Conduction.Functions.Dw(
      w[i - 1],
      material.w_f,
      material.A_layer,
      material.dww_tab,
      material.dws_tab,
      switch_dw,
      activatesuction)));
  end for;
  D_phi[nSta + 1] = dw_dphi[nSta]*
    Buildings.HeatTransfer.Conduction.Functions.Dw(
    w[nSta],
    material.w_f,
    material.A_layer,
    material.dww_tab,
    material.dws_tab,
    switch_dw,
    activatesuction)/(d_node/2);

  //THERMAL CONDUCTANCE
  UA[1] = A*lambda[1]/(d_node/2);
  for i in 2:nSta loop
    UA[i] = 1.0/((d_node/2)/(lambda[i]*A) + (d_node/2)/(lambda[i - 1]*A));
  end for;
  UA[nSta + 1] = lambda[nSta]/(d_node/2);

  //TRANSPORT EQUATIONS
  for i in 2:nSta loop
    T[i - 1] - T[i] = Q_flow[i]/UA[i];
    (pw[i - 1] - pw[i])*(kd[i]*A) + (phi[i - 1] - phi[i])*D_phi[i] = m_flow[i];

  end for;

  for i in 1:nSta loop
    Cm[i] = Cd + w[i]*Cw*V;
    psat[i] = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat(T[i]);
    phi[i] = pw[i]/psat[i];
    Cm[i]*der(T[i]) + V*Cw*T[i]*dw_dphi[i]*der(phi[i]) = (Q_flow[i] - Q_flow[i + 1]) + hv*(m_flow[i] - m_flow[i + 1]);
    dw_dphi[i]*der(phi[i]) = (m_flow[i] - m_flow[i + 1])/V;
  end for;

  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-96,4},{96,-4}},
          lineColor={0,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.HorizontalCylinder),
        Polygon(
          points={{12,14},{14,14},{16,10},{18,4},{18,0},{14,-6},{12,-8},{6,-14},
              {-2,-16},{-6,-12},{-12,-6},{-14,4},{-12,10},{-10,14},{-8,16},{-4,
              18},{-2,18},{4,18},{8,18},{12,14}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-6,-10},{2,-14},{8,-10},{14,-6},{10,-10},{8,-16},{-2,-18},{-10,
              -14},{-12,-6},{-14,4},{-12,10},{-10,14},{-8,16},{-10,6},{-10,-2},
              {-6,-10}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,72},{-44,-72}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{44,72},{60,-72}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Line(
          points={{-44,40},{44,40}},
          color={0,0,255},
          pattern=LinePattern.DashDotDot,
          smooth=Smooth.None,
          arrow={Arrow.Filled,Arrow.Filled},
          thickness=0.5),
        Line(
          points={{-44,-40},{44,-40}},
          color={0,0,255},
          pattern=LinePattern.DashDotDot,
          smooth=Smooth.None,
          arrow={Arrow.Filled,Arrow.Filled},
          thickness=0.5),
        Text(
          extent={{-106,-70},{0,-88}},
          lineColor={0,0,255},
          textString="%x"),
        Text(
          extent={{-36,-64},{42,-94}},
          lineColor={0,0,255},
          textString="%nSta")}),
    defaultComponentName="lay",
    Documentation(info="<html>
    This is a model of a heat and mass conductor for a single layer of homogeneous material
    that computes transient  heat and moisture conduction. The material is a record that extends
<a href=\"modelica://Buildings.HeatTransfer.Data.Solids\">
Buildings.HeatTransfer.Data.Solids.</a> 
    

<h4>Transient heat conduction in materials </h4>

</p>
<p align=\"center\" style=\"font-style:italic;\">
   (dH &frasl; dT)&sdot;(&part; T(s,t) &frasl; &part;t) = 
   &part;(&lambda;(w)&sdot; &part;T(s,t) &frasl; &part;s)/&part;s + h<sub>v</sub>&sdot;&part;( &delta;<sub>p</sub>&sdot; (&part;(&phi;&sdot;P<sub>sat</sub>)&frasl; &part;s)) &frasl; &part;s
</p>
<p>
Where : <ul>
<p> <i>dH &frasl; dT</i>
is the heat storage capacity of the moist building material in [J/m&sup3;K]</P>
<p>H is the enthalpy of the building material in [J/m<sup>3</sup>] : 
<p align=\"center\" style=\"font-style:italic;\">
H = &rho;c<sub>s</sub>&sdot;T + c<sub>w</sub>&sdot;w
</p>
</p>
<p>
<i>T</i>
is the temperature in [K] at location <i>s</i> and time <i>t</i> </p>
<p>
<i>h<sub>v</sub></i>   is the evaporation enthalpy of the water in [J/kg]</p> 
<p> <i>P<sub>sat</sub></i>  the water vapour saturation pressure [Pa]
</p>
<p>
&phi; is the relative humidity 
</p>
<p>
<i>&lambda;(w)</i>
is the thermal conductivity of moist building material in [W/mK], it can be determined by two different ways :
</p>

<ul><p >
 If <b> switch_lamb = 1 </b>then <i>&lambda; = &lambda;&sdot;(1 + b_h&sdot;w/&rho;)</i>
</p>
<p>
<i>b_h</i>
is the thermal conductivity supplement, an approximation factor stored in the material hygro thermal properties and used to determine the influence of moisture on the thermal conductivity. 
</p>
<p>
 If <b> switch_lamb is equal to any other integer </b> then  &lambda; is determined by a linear interpolation of the table stored in the material record using 
the fonction <a href=\"modelica://Modelica.Math.Vectors.interpolate\"> Modelica.Math.Vectors.interpolate </a>
</p></ul>
<p>
<i>w</i>  is the water content of the building material in [kg/m<sup>3</sup>] at location <i>s</i> and time <i>t</i>, it depend of the sorption isotherm of the material :
</p>
<ul>  If <b> switch_w = 1 </b> then  w is a determined by an approximation of the sorption isotherm :
<p align=\"center\" style=\"font-style:italic;\">w = w<sub>f</sub>&sdot;((b - 1)&sdot;&phi;/(b - &phi;)) </p>
<p align=\"center\" style=\"font-style:italic;\">b = 0.8&sdot;(w<sub>80</sub>-w<sub>f</sub>)/(w<sub>80</sub>-0.8&sdot;w<sub>f</sub>)</p>
If <b> switch_w = 2 </b> then w is determined interpolating the table of the sorption isotherm stored in the material record.
<p> If <b> switch_w is equal to any other integer </b> w is determined using an approximation of the sorption isotherm that does not require specific hygrothermal properties of the material
but only the porosity of this particular material.</p> </ul>


   <p>
   <i>&delta;<sub>p</sub></i> is the water permeability of the building material [kg/msPa] :
<p align=\"center\" style=\"font-style:italic;\">
&delta;(T) = 2.0&sdot;10<sup>-7</sup>T<sup>0.81</sup> &frasl; (P<sub>L</sub>&sdot;&mu;)
</p>
<p>
Remark : h<sub>v</sub>&sdot;&part;( &delta;<sub>p</sub>&sdot; (&part;(&phi;&sdot;P<sub>sat</sub>)&frasl; &part;s)) &frasl; &part;s is the influence of water on the heat transfert.
</p>


<h4>Transient moisture conduction in materials </h4>

</p>
<p align=\"center\" style=\"font-style:italic;\">
   (dw &frasl; d&phi;)(&part; &phi; &frasl; &part;t) = 
   &part; (D&#7529;&sdot;(&part; &phi; &frasl; &part;s) + &delta;<sub>p</sub>&sdot; (&part; (&phi;&sdot;P<sub>sat</sub>)&frasl; &part;s) &frasl; &part;s
</p>
<p>
Where 
<i>dw &frasl; d&phi;</i> is the moisture storage capacity of the building material in [kg/m<sup>3</sup>], as for the water content w, it is determined by three diffenrent ways using the same switch
(switch = 1 Kunzel's approximation, switch = 2 Interpolation of the derivative of the sorption isotherm, switch = any other integer: approximation without known sorption isotherm).
</p>

<p>
D&#7529; is the liquid conduction coefficient in [kg/ms] :
<p align=\"center\" style=\"font-style:italic;\">
D&#7529; = D<sub>w</sub>&sdot;(dw &frasl; d&phi;)
</p>
<p>
Dw is the capillarity transport coefficient [m<sup>2</sup>], if the material is in suction process (parameter <var>activatesuction</var> should be use when the building component is in contact with water for short periods of time, e.g. when it rains) D<sub>w</sub> = D<sub>ws</sub> if not D<sub>w</sub> = D<sub>ww</sub>.
<ul> If <b> switch_Dw = 1 </b> <p align=\"center\" style=\"font-style:italic;\"> D<sub>ws</sub> = 3.8&sdot;(A/w<sub>f</sub>)<sup>2</sup>&sdot;1000<sup>w/w<sub>f</sub>-1</sup> </p>
<p align=\"center\" style=\"font-style:italic;\"> D<sub>ww</sub> = D<sub>ws</sub>/10 </p>

If <b> switch_Dw = 2 </b> D<sub>ww</sub> and D<sub>ws</sub> are interpolated from the tables stored in the material record. </ul>


<h4>Spatial discretization</h4>
<p>
To spatially discretize the heat and moisture equations, the construction is 
divided into compartments with <code>material.nSta &ge; 1</code> state variables. 
The state variables are connected to each other through hygrothermal conductors. 
There is also an hygrothermal conductor
between the surfaces and the outermost state variables. Thus, to obtain
the surface temperature, use <code>heatMassPort_a.heatPort.T</code> (or <code>heatMassPort_b.heatPort.T</code>)
and not the variable <code>T[1]</code></p>
<h4>References</h4>


<p>
Hartwig M. Kunzel. 1995. <u>Simultaneous Heat and Moisture transport in Building Components.</u> Franhofer Institute of Building Physics (ISBN 3-8167-4103-7).
</p>
</html>", revisions="<html>
<ul>
<li>
Jun 8 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
end SingleLayerHM;
