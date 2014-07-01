within Buildings.HeatTransfer.Convection.BaseClasses;
partial model PartialConvectionHM
  "Partial model for heat and water vapour convection"
import Buildings;
extends Buildings.BaseClasses.BaseIcon;
parameter Modelica.SIunits.Area A "Heat transfer area";
parameter Modelica.SIunits.CoefficientOfHeatTransfer hFixed=3
    "Constant convection coefficient";

Modelica.SIunits.HeatFlowRate Q_flow "Heat flow rate from solid -> fluid";
Modelica.SIunits.HeatFlux q_flow "Convective heat flux from solid -> fluid";
Modelica.SIunits.TemperatureDifference dT(start=0) "= solid.T - fluid.T";
Modelica.SIunits.MassFlowRate m_flow "Mass flow rate from solid -> fluid";
Real g_flow "Water vapour flux density [kg.m^-2.s^-1]";
Modelica.SIunits.PressureDifference dP;
Modelica.SIunits.Pressure P_fluid;
Modelica.SIunits.Pressure P_solid;

parameter Modelica.SIunits.Angle til(displayUnit="deg") "Surface tilt"
  annotation (Dialog(enable= not (conMod == Buildings.HeatTransfer.Types.InteriorConvection.fixed)));
protected
final parameter Real cosTil=Modelica.Math.cos(til) "Cosine of window tilt"
  annotation (Evaluate=true);
final parameter Real sinTil=Modelica.Math.sin(til) "Sine of window tilt"
  annotation (Evaluate=true);
final parameter Boolean isCeiling = abs(sinTil) < 10E-10 and cosTil > 0
    "Flag, true if the surface is a ceiling"
  annotation (Evaluate=true);
final parameter Boolean isFloor = abs(sinTil) < 10E-10 and cosTil < 0
    "Flag, true if the surface is a floor"
  annotation (Evaluate=true);

public
Buildings.HeatTransfer.Interfaces.HeatMassPort_a solid annotation (
    Placement(transformation(extent={{-110,-4},{-90,16}}), iconTransformation(
        extent={{-110,-4},{-90,16}})));
Buildings.HeatTransfer.Interfaces.HeatMassPort_b fluid annotation (
    Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(
        extent={{90,-10},{110,10}})));
equation
dT = solid.heatPort.T - fluid.heatPort.T;
solid.heatPort.Q_flow = Q_flow;
fluid.heatPort.Q_flow = -Q_flow;
Q_flow = A*q_flow;
P_solid = Buildings.Utilities.Psychrometrics.Functions.pW_X(solid.massPort.Xi_outflow);
P_fluid = Buildings.Utilities.Psychrometrics.Functions.pW_X(fluid.massPort.Xi_outflow);

dP = P_solid-P_fluid;

solid.massPort.m_flow = m_flow;
fluid.massPort.m_flow = -m_flow;
m_flow = A*g_flow;

annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
          -100},{100,100}}), graphics={
      Rectangle(
        extent={{-100,100},{100,-100}},
        lineColor={0,0,0},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid),
      Rectangle(
        extent={{-80,90},{-50,-70}},
        lineColor={0,0,0},
        fillColor={192,192,192},
        fillPattern=FillPattern.Backward),
      Line(points={{-24,90},{-24,-70}}, color={0,127,255}),
      Line(points={{16,90},{16,-70}},
                                    color={0,127,255}),
      Line(points={{50,90},{50,-70}}, color={0,127,255}),
      Line(points={{86,90},{86,-70}}, color={0,127,255}),
      Line(points={{-24,-70},{-14,-50}}, color={0,127,255}),
      Line(points={{16,-70},{26,-50}},   color={0,127,255}),
      Line(points={{50,-70},{60,-50}},   color={0,127,255}),
      Line(points={{86,-70},{96,-50}},   color={0,127,255}),
      Line(points={{-24,-70},{-34,-50}}, color={0,127,255}),
      Line(points={{16,-70},{6,-50}},    color={0,127,255}),
      Line(points={{50,-70},{40,-50}},   color={0,127,255}),
      Line(points={{86,-70},{76,-50}},   color={0,127,255}),
      Line(points={{-50,30},{86,30}}, color={191,0,0}),
      Line(points={{66,40},{86,30}}, color={191,0,0}),
      Line(points={{66,20},{86,30}}, color={191,0,0}),
      Text(
        extent={{-25,52},{5,30}},
        lineColor={255,0,0},
        textString="Q_flow"),
      Line(points={{-50,-10},{86,-10}}, color={0,127,0}),
      Line(points={{66,0},{86,-10}},   color={0,127,0}),
      Line(points={{66,-20},{86,-10}}, color={0,127,0}),
      Text(
        extent={{-25,-10},{5,-32}},
        lineColor={0,128,0},
        textString="m_flow")}));
end PartialConvectionHM;
