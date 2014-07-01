within Buildings.HeatTransfer.Conduction;
model SingleLayerHM

  Modelica.SIunits.MassFlowRate m_flow[nSta + 1];
  Modelica.SIunits.HeatFlowRate Q_flow[nSta + 1];
  Modelica.SIunits.Temperature T[nSta]( each start=T_ini)
    "Temperature at the states";

  Real kd[nSta + 1];

  Modelica.SIunits.Pressure pw[nSta] "Vapour Pressure at the state";

  Real phi[nSta](each start=phi_ini) "Relative humidity";

  Real dw_dphi[nSta] "Derivative of water content over relative humidity.";

  Modelica.SIunits.Pressure psat[nSta]
    "equilibrum vapor pressure at the states";

  Modelica.SIunits.MassFraction Xi_outflow[nSta] "mass fraction at the state";

  Real w[nSta] "Water content at the states";

  Modelica.SIunits.HeatCapacity Cm[nSta]
    "Heat capacity of moist material at the states";

  Modelica.SIunits.ThermalConductivity lambda[nSta];

  Modelica.SIunits.ThermalConductance UA[nSta + 1];

  Real D_phi[nSta + 1] "liquid conduction coefficient [kg/ms]";

    parameter Modelica.SIunits.Area A=1 "Heat transfer area";
    parameter Real phi_ini = 0.0
    "Initial value of relative humidity at the nodes";
    parameter Modelica.SIunits.Temperature T_ini = 293.15;

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
  parameter Real b=0.8*(material.w_80-material.w_f)/(material.w_80-0.8*material.w_f);

  parameter Modelica.SIunits.Temperature T_a_start=293.15
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  parameter Modelica.SIunits.Temperature T_b_start=293.15
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  constant Modelica.SIunits.Pressure patm=101325 "Atmospheric pressure";
  constant Modelica.SIunits.SpecificHeatCapacity Cw=4185.0;
  constant Real hv=2500000 "Evaporation heat of water [J/Kg]";

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
    = m_flow[nSta + 1]/(kd[nSta+1]*A);

  for i in 1:nSta loop
    dw_dphi[i] =Buildings.HeatTransfer.Conduction.Functions.dw_dphi(
      phi[i],
      material.w_f,
      material.w_80,
      b,
      material.sorp_tab_layer,
      material.Kunzel);

    lambda[i]  =Buildings.HeatTransfer.Conduction.Functions.lambdaHM(
      material.d,
      material.k,
      w[i],
      b,
      material.lamb_tab_layer,
      material.Kunzel);

    w[i]       =Buildings.HeatTransfer.Conduction.Functions.sorption(
      phi[i],
      material.w_f,
      material.w_80,
      b,
      material.sorp_tab_layer,
      material.Kunzel);

    Xi_outflow[i] = Buildings.Utilities.Psychrometrics.Functions.pW_X(pw[i]);

 end for;

  // WATER VAPOUR DIFFUSION COEFFICIENT
  kd[1] = Buildings.HeatTransfer.Conduction.Functions.delta_L(T[1])/material.mu/(material.x/(2 *nSta));
  for i in 2:nSta loop
     kd[i] = 1.0 / ((material.x/(2 *nSta))/(Buildings.HeatTransfer.Conduction.Functions.delta_L(T[i])/material.mu)
                  + (material.x/(2 *nSta))/(Buildings.HeatTransfer.Conduction.Functions.delta_L(T[i-1])/material.mu));
  end for;
  kd[nSta+1] = Buildings.HeatTransfer.Conduction.Functions.delta_L(T[nSta])/material.mu/(material.x/(2*nSta));

  // LIQUID WATER DIFFUSION COEFFIIENT
  D_phi[1] = dw_dphi[1] * (Buildings.HeatTransfer.Conduction.Functions.Capillary_transp_coef(w[1],
                                                                                             material.w_f,
                                                                                             material.A,
                                                                                             material.dww_tab_layer,
                                                                                             material.Kunzel))/(material.x/(2 *nSta));
  for i in 2:nSta loop
    D_phi[i] =  1.0 / ((material.x/(2 *nSta)) / (Buildings.HeatTransfer.Conduction.Functions.Capillary_transp_coef(w[i],
                                                                                                    material.w_f,
                                                                                                    material.A,
                                                                                                    material.dws_tab_layer,
                                                                                                    material.Kunzel)))
                        + ((material.x/(2 *nSta)) / (dw_dphi[1] * (Buildings.HeatTransfer.Conduction.Functions.Capillary_transp_coef(w[i],
                                                                                                    material.w_f,
                                                                                                    material.A,
                                                                                                    material.dww_tab_layer,
                                                                                                    material.Kunzel))));
  end for;
  D_phi[nSta+1] = dw_dphi[nSta] * (Buildings.HeatTransfer.Conduction.Functions.Capillary_transp_coef(w[nSta],
                                                                                                    material.w_f,
                                                                                                    material.A,
                                                                                                    material.dww_tab_layer,
                                                                                                    material.Kunzel))/(material.x/(2 *nSta));

  //THERMAL CONDUCTANCE
  UA[1] = lambda[1]*A/(material.x/(2 *nSta));
  for i in 2:nSta loop
    UA[i] = 1.0 /  ( (material.x/(2 *nSta))/(lambda[i]*A)
                  +  (material.x/(2 *nSta))/(lambda[i - 1]*A));
  end for;
  UA[nSta+1] = lambda[nSta]*A/(material.x/(2*nSta));

  //TRANSPORT EQUATIONS
  for i in 2:nSta loop
    T[i - 1] - T[i] = Q_flow[i]/UA[i];
    (pw[i - 1] - pw[i])* (kd[i]*A) +  (phi[i - 1] - phi[i]) * D_phi[i] = m_flow[i];

  end for;

  for i in 1:nSta loop
    Cm[i] = Cd + w[i]*Cw*V;
    psat[i]= Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat(T[i]);
    phi[i] = pw[i]/psat[i];
    Cm[i]*der(T[i]) + V*Cw*T[i]*dw_dphi[i]*der(phi[i]) = (Q_flow[i] - Q_flow[i + 1])+hv*(m_flow[i] - m_flow[i + 1])/V;
    dw_dphi[i]*der(phi[i]) = (m_flow[i] - m_flow[i + 1])/V;
  end for;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(coordinateSystem(preserveAspectRatio=false, extent=
            {{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-96,4},{96,-4}},
          lineColor={0,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.HorizontalCylinder),
        Polygon(
          points={{12,14},{14,14},{16,10},{18,4},{18,0},{14,-6},{12,-8},{6,-14},{-2,-16},{-6,-12},
              {-12,-6},{-14,4},{-12,10},{-10,14},{-8,16},{-4,18},{-2,18},{4,18},{8,18},{12,14}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-6,-10},{2,-14},{8,-10},{14,-6},{10,-10},{8,-16},{-2,-18},{-10,-14},{-12,-6},{
              -14,4},{-12,10},{-10,14},{-8,16},{-10,6},{-10,-2},{-6,-10}},
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
    that computes transient  heat and moisture conduction. If the material is a record that extends
<a href=\"modelica://Buildings.HeatTransfer.Data.Solids\">
Buildings.HeatTransfer.Data.Solids</a> and its parameter Kunzel is true then the transient combined heat and moisture transfer is calculated using the parameters described in <u> Simultaneous Heat and Moisture Transport in Building Components </u> Hartwig M. Kunzel (Fraunhofer Institute of Building Physics) else the model uses the parameters stored in the material record.
    

<h4>Transient heat conduction in materials </h4>

</p>
<p align=\"center\" style=\"font-style:italic;\">
   (dH &frasl; dT)&sdot;(&part; T(s,t) &frasl; &part;t) = 
   &part;(&lambda;(w)&sdot; &part;T(s,t) &frasl; &part;s)/&part;s + h<sub>v</sub>&sdot;&part;( &delta;<sub>p</sub>&sdot; (&part;(&phi;&sdot;P<sub>sat</sub>)&frasl; &part;s)) &frasl; &part;s
</p>
<p>
where 
<i>dH &frasl; dT</i>
is the heat storage capacity of the moist building material in [J/m&sup3;K], H is the enthalpy of the building material in [J/m<sup>3</sup>] : 
<p align=\"center\" style=\"font-style:italic;\">
H = &rho;c<sub>s</sub>&sdot;T + c<sub>w</sub>&sdot;w
   </p>
<p>
<i>T</i>
is the temperature in [K] at location <i>s</i> and time <i>t</i> , <i>w</i> is the water content of the building material at location <i>s</i> and time <i>t</i> in [kg/m<sup>3</sup>], <i>h<sub>v</sub></i>  is the evaporation enthalpie of the water and <i>P<sub>sat</sub></i> the water vapour saturation pressure.
</p>
<i>&lambda;(w)</i>
is the thermal conductivity of moist building material in [W/mK]
</p>

<p align=\"center\" style=\"font-style:italic;\">
&lambda;(w) = &lambda;(1 + b&sdot;w/&rho;)
</p>
<p>
<i>b</i>
[.] is an approximation factor used by Kunzel to take into account the influence of moisture on the thermal conductivity of the material. It must be greater than one and it can be determined from the equilibrium water content at 80% of relative humidity.
<p align=\"center\" style=\"font-style:italic;\">
b = 0.8 (w_80-w_f)/(w_80-0.8*w_f)
</p>
If Kunzel is false, the parameter b is not used in the model, instead the model uses linear interpolations of the sorption isotherm with the data stored in <a href=\"modelica://Buildings.HeatTransfer.Data.BaseClasses.HygroThermalMaterial\"> Buildings.HeatTransfer.Data.BaseClasses.HygroThermalMaterial </a>. If Kunzel is true, b is used to calculate dw &frasl; d&phi;, w and &lambda; at each node.
   <p>
   <i>&delta;<sub>p</sub></i> is the water permeability of the building material [kg/msPa]
   <p align=\"center\" style=\"font-style:italic;\">
&delta;<sub>p</sub> = &delta;/&mu;
</p>

 <p>  
   <i>&delta;</i> is the water vapour diffusion coefficient in air [kg/msPa]
<p align=\"center\" style=\"font-style:italic;\">
&delta;(T) = 2.0&sdot;10<sup>-7</sup>T<sup>0.81</sup> &frasl; P<sub>L</sub>
</p>
<p>
h<sub>v</sub>&sdot;&part;( &delta;<sub>p</sub>&sdot; (&part;(&phi;&sdot;P<sub>sat</sub>)&frasl; &part;s)) &frasl; &part;s is the influence of water on the heat transfert.
</p>


<h4>Transient moisture conduction in materials </h4>

</p>
<p align=\"center\" style=\"font-style:italic;\">
   (dw &frasl; d&phi;)(&part; &phi; &frasl; &part;t) = 
   &part; (D&#7529;&sdot;(&part; &phi; &frasl; &part;s) + &delta;<sub>p</sub>&sdot; (&part; (&phi;&sdot;P<sub>sat</sub>)&frasl; &part;s) &frasl; &part;s
</p>
<p>
where 
<i>dw &frasl; d&phi;</i> is the moisture storage capacity of the building material in [kg/m<sup>3</sup>],
<i>&phi;</i> is the relative humidity,

<p>
D&#7529; is the liquid conduction coefficient in [kg/ms]
<p align=\"center\" style=\"font-style:italic;\">
D&#7529; = D<sub>w</sub>&sdot;(dw &frasl; d&phi;)
</p>
Idem that for the heat conduction Dw can be interpolated from the data stored in the record or calculated using Kunzel's approximation. 
<h4>Spatial discretization</h4>
<p>
To spatially discretize the heat and moisture equations, the construction is 
divided into compartments with <code>material.nSta &ge; 1</code> state variables. 
The state variables are connected to each other through hygrothermal conductors. 
There is also an hygrothermal conductor
between the surfaces and the outermost state variables. Thus, to obtain
the surface temperature, use <code>heatMassPort_a.heatPort.T</code> (or <code>heatMassPort_b.heatPort.T</code>)
and not the variable <code>T[1]</code></p>

          </html>"));
end SingleLayerHM;
