within Buildings.Fluid.HeatExchangers.DXCoils.BaseClasses.Examples;
model DryCoil "Test model for DryCoil"
extends Modelica.Icons.Example;
  package Medium =
      Buildings.Media.GasesConstantDensity.MoistAirUnsaturated;
  Modelica.Blocks.Sources.Constant p(
    k=101325)
    annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
  Modelica.Blocks.Sources.Constant TConIn(
    k=273.15 + 35) "Condenser inlet air temperature"
    annotation (Placement(transformation(extent={{-80,46},{-60,66}})));
  Buildings.Fluid.HeatExchangers.DXCoils.BaseClasses.DryCoil dryCoi(
    redeclare package Medium = Medium,
    datCoi=datCoi) "Performs calculation for dry coil condition"
    annotation (Placement(transformation(extent={{20,-10},{40,10}})));
  Modelica.Blocks.Sources.BooleanStep onOff(
    startTime=1200) "Compressor on-off signal"
    annotation (Placement(transformation(extent={{-20,20},{0,40}})));
  Modelica.Blocks.Sources.Ramp m_flow(
    startTime=1200,
    duration=600,
    height=1.5) "Mass flow rate of air"
    annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
  Modelica.Blocks.Sources.Ramp TIn(
    duration=600,
    startTime=2400,
    height=-4,
    offset=273.15 + 29) "Dry bulb temperature of air entring the coil"
    annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
  Modelica.Blocks.Sources.Ramp XIn(
    duration=600,
    startTime=2400,
    height=-0.002,
    offset=0.006) "Inlet mass-fraction"
    annotation (Placement(transformation(extent={{-80,-94},{-60,-74}})));
  Modelica.Blocks.Sources.Ramp hIn(
    duration=600,
    startTime=2400,
    height=-10000,
    offset=45000) "Specific enthalpy of air entring the coil"
    annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
  Modelica.Blocks.Sources.TimeTable speRat(table=[0.0,0.0; 900,0.25; 1800,0.50;
        2700,0.75]) "Speed ratio "
    annotation (Placement(transformation(extent={{-80,76},{-60,96}})));
  Data.CoilData datCoi(per={
        Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.Generic(
        spe=900,
        nomVal=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.NominalValues(
          Q_flow_nominal=-12000,
          COP_nominal=3,
          SHR_nominal=0.8,
          m_flow_nominal=0.9),
        perCur=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.PerformanceCurves.Curve_I()),
        Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.Generic(
        spe=1200,
        nomVal=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.NominalValues(
          Q_flow_nominal=-18000,
          COP_nominal=3,
          SHR_nominal=0.8,
          m_flow_nominal=1.2),
        perCur=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.PerformanceCurves.Curve_I()),
        Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.Generic(
        spe=1800,
        nomVal=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.NominalValues(
          Q_flow_nominal=-21000,
          COP_nominal=3,
          SHR_nominal=0.8,
          m_flow_nominal=1.5),
        perCur=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.PerformanceCurves.Curve_II()),
        Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.Generic(
        spe=2400,
        nomVal=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.BaseClasses.NominalValues(
          Q_flow_nominal=-30000,
          COP_nominal=3,
          SHR_nominal=0.8,
          m_flow_nominal=1.8),
        perCur=
          Buildings.Fluid.HeatExchangers.DXCoils.Data.PerformanceCurves.Curve_III())}, nSpe
      =4) "Coil data"
    annotation (Placement(transformation(extent={{60,60},{80,80}})));
equation
  connect(TConIn.y, dryCoi.TConIn)  annotation (Line(
      points={{-59,56},{-42,56},{-42,5},{19,5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(p.y, dryCoi.p)  annotation (Line(
      points={{-59,-50},{-42,-50},{-42,-2.4},{19,-2.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(onOff.y, dryCoi.on)  annotation (Line(
      points={{1,30},{12,30},{12,10},{19,10}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(m_flow.y, dryCoi.m_flow)  annotation (Line(
      points={{-59,22},{-48,22},{-48,2.4},{19,2.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TIn.y, dryCoi.TIn)  annotation (Line(
      points={{-59,-18},{-48,-18},{-48,6.10623e-16},{19,6.10623e-16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(XIn.y, dryCoi.XIn)  annotation (Line(
      points={{-59,-84},{-36,-84},{-36,-5},{19,-5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(hIn.y, dryCoi.hIn)  annotation (Line(
      points={{1,-30},{10,-30},{10,-7.7},{19,-7.7}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(speRat.y, dryCoi.speRat)     annotation (Line(
      points={{-59,86},{-36,86},{-36,7.6},{19,7.6}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Fluid/HeatExchangers/DXCoils/BaseClasses/Examples/DryCoil.mos"
        "Simulate and plot"),
          Documentation(info="<html>
<p>
This example illustrates working of DryCoil block 
<a href=\"modelica://Buildings.Fluid.HeatExchangers.DXCoils.BaseClasses.DryCoil\">
Buildings.Fluid.HeatExchangers.DXCoils.BaseClasses.DryCoil</a>. 
</p>
</html>",
revisions="<html>
<ul>
<li>
April 10, 2012 by Kaustubh Phalak:<br>
First implementation. 
</li>
</ul>
</html>"));
end DryCoil;