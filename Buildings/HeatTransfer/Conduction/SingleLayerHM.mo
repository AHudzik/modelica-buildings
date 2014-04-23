within Buildings.HeatTransfer.Conduction;
model SingleLayerHM

  Modelica.SIunits.MassFlowRate m_flow[nSta + 1];
  Modelica.SIunits.HeatFlowRate Q_flow[nSta + 1];
  Modelica.SIunits.Temperature T[nSta]( each start=T_ini)
    "Temperature at the states";
  parameter Modelica.SIunits.Temperature T_ini = 293.15;
    Real kd[nSta + 1];
  //Real D[nSta]
  //  "Water vapour diffusion coefficient in air (Unit = Kg.m^-2.s^-1)";
  Modelica.SIunits.Pressure pw[nSta] "Vapour Pressure at the state";
  Real phi[nSta](each start=phi_ini) "Relative humidity";
  parameter Real phi_ini = 0.5
    "initial value of relative humidity at the nodes";
  Real dw_dphi[nSta] "Derivative of water content over relative humidity.";

  Modelica.SIunits.Pressure psat[nSta]
    "equilibrum vapor pressure at the states";

  //Modelica.SIunits.MassFraction Xi_outflow[nSta] "mass fraction at the state";
  Modelica.SIunits.MassConcentration w[nSta] "Water content at the states";
  Modelica.SIunits.HeatCapacity Cm[nSta]
    "Heat capacity of moist material at the states";
    Modelica.SIunits.ThermalConductivity lambda[nSta];
    Modelica.SIunits.ThermalConductance UA[nSta + 1];
    Real D_phi[nSta + 1] "liquid conduction coefficient [kg/ms]";

  replaceable parameter Data.BaseClasses.HygroThermalMaterial material
    "Material from Data.Solids, Data.SolidsPCM or Data.Resistances" annotation (
    Evaluate=true,
    choicesAllMatching=true,
    Placement(transformation(extent={{60,60},{80,80}})));

protected
  parameter Modelica.SIunits.Area A=1 "Heat transfer area";
  //final parameter Modelica.SIunits.CoefficientOfHeatTransfer U=UA/A
    //"U-value (without surface heat transfer coefficients)";
  //parameter Modelica.SIunits.ThermalResistance R=material.x/(material.k*A)
   // "Thermal resistance of construction";
  //final parameter Modelica.SIunits.ThermalConductance UA=1/R
   // "Thermal conductance of construction (without surface heat transfer coefficients)";

  Modelica.SIunits.TemperatureDifference dT "port_a.T - port_b.T";
  //final parameter Real l=material.x/nSta;
  final parameter Modelica.SIunits.Volume V=A*material.x/nSta;
  final parameter Integer nSta(min=1) = material.nSta
    "Number of state variables";
 // final parameter Modelica.SIunits.ThermalConductance UAnSta=UA*nSta
   // "Thermal conductance between nodes";
 // final parameter Modelica.SIunits.ThermalConductance UAnSta2=2*UAnSta
   // "Thermal conductance between nodes and surface boundary";
  parameter Modelica.SIunits.Mass m=A*material.x*material.d/material.nSta
    "Mass associated with the temperature state";
  parameter Modelica.SIunits.HeatCapacity Cd=m*material.c
    "Heat capacity of dry material associated with the temperature state";
  parameter Real b=material.b;

  parameter Modelica.SIunits.Temperature T_a_start=293.15
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  parameter Modelica.SIunits.Temperature T_b_start=293.15
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  constant Modelica.SIunits.Pressure patm=101325 "Atmospheric pressure";
  constant Modelica.SIunits.SpecificHeatCapacity Cw=4185.0;

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
    dw_dphi[i] = material.w_f*(b - 1)*b/((b - phi[i])^2);
    w[i] = material.w_f*((b - 1)*phi[i])/(b - phi[i]);
    lambda[i]=material.k*(1+b*w[i]/material.d);
  end for;

  kd[1] = Buildings.HeatTransfer.Conduction.Functions.delta_L(T[1])/material.mu/(material.x/(2 *nSta));
  for i in 2:nSta loop
     kd[i] = 1.0 / ((material.x/(2 *nSta))/(Buildings.HeatTransfer.Conduction.Functions.delta_L(T[i])/material.mu)
                  + (material.x/(2 *nSta))/(Buildings.HeatTransfer.Conduction.Functions.delta_L(T[i-1])/material.mu));
  end for;
  kd[nSta+1] = Buildings.HeatTransfer.Conduction.Functions.delta_L(T[nSta])/material.mu/(material.x/(2*nSta));

  D_phi[1] = dw_dphi[1] * (3.8*(material.A/material.w_f)^2*1000^(w[1]/(material.w_f - 1)))/(material.x/(2 *nSta));
  for i in 2:nSta loop
    D_phi[i] =  1.0 / ((material.x/(2 *nSta)) / (dw_dphi[1] * (3.8*(material.A/material.w_f)^2*1000^(w[1]/(material.w_f - 1))))
                        + ((material.x/(2 *nSta)) / (dw_dphi[1] * (3.8*(material.A/material.w_f)^2*1000^(w[1]/(material.w_f - 1))))));
  end for;

  D_phi[nSta+1] = dw_dphi[nSta] * (3.8*(material.A/material.w_f)^2*1000^(w[nSta]/(material.w_f - 1)))/(material.x/(2 *nSta));

  UA[1] = lambda[1]*A/(material.x/(2 *nSta));
  for i in 2:nSta loop
    UA[i] = 1.0 /  ( (material.x/(2 *nSta))/(lambda[i]*A)
                  +  (material.x/(2 *nSta))/(lambda[i - 1]*A));
  end for;
  UA[nSta+1] = lambda[nSta]*A/(material.x/(2*nSta));

  for i in 2:nSta loop
    T[i - 1] - T[i] = Q_flow[i]/UA[i];
    (pw[i - 1] - pw[i])* (kd[i]*A) +  (phi[i - 1] - phi[i]) * D_phi[i] = m_flow[i];

  end for;

  for i in 1:nSta loop
    Cm[i] = Cd + w[i]*Cw*V;
    psat[i]= Buildings.HeatTransfer.Conduction.Functions.p_sat(T[i]);
    phi[i] = pw[i]/psat[i];
    Cm[i]*der(T[i]) + V*Cw*T[i]*dw_dphi[i]*der(phi[i]) = (Q_flow[i] - Q_flow[i + 1])+2500000*(m_flow[i] - m_flow[i + 1])/V;
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
          textString="%nSta")}));

end SingleLayerHM;
