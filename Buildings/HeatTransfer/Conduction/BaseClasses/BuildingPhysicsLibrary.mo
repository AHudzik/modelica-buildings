within Buildings.HeatTransfer.Conduction.BaseClasses;
model BuildingPhysicsLibrary
  extends Modelica.Icons.Example;
  SingleLayerHM lay(
    phi_ini=0,
    A=1,
    switch_lamb=1,
    T_ini(displayUnit="degC") = 278.15,
    switch_w=1,
    switch_dw=1,
    activatesuction=true,
    redeclare Buildings.HeatTransfer.Data.Solids.GypsumBoardHM material(x=0.1,
        nSta=8))
    annotation (Placement(transformation(extent={{-18,-18},{16,18}})));
  Interfaces.HeatMassPort_b heatMassPort_b1
    annotation (Placement(transformation(extent={{36,-4},{44,4}})));
  Interfaces.HeatMassPort_a heatMassPort_a1
    annotation (Placement(transformation(extent={{-44,-4},{-36,4}})));
  Sources.FixedHumidity fixedHumidity(X=0.0015)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Sources.FixedTemperature fixedTemperature(T(displayUnit="degC") = 278.15)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Sources.FixedTemperature fixedTemperature1(T=303.15)
    annotation (Placement(transformation(extent={{82,20},{62,40}})));
  Sources.FixedHumidity fixedHumidity1(X=0.0218)
    annotation (Placement(transformation(extent={{80,-40},{60,-20}})));
equation
  connect(lay.heatMassPort_b, heatMassPort_b1) annotation (Line(
      points={{16,0},{40,0}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(lay.heatMassPort_a, heatMassPort_a1) annotation (Line(
      points={{-18,0},{-40,0}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(fixedTemperature.port, heatMassPort_a1.heatPort) annotation (Line(
      points={{-60,30},{-50,30},{-50,0},{-40,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
      points={{-60,-30},{-50,-30},{-50,0},{-40,0}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature1.port, heatMassPort_b1.heatPort) annotation (Line(
      points={{62,30},{52,30},{52,0},{40,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity1.massPort, heatMassPort_b1.massPort) annotation (Line(
      points={{60,-30},{52,-30},{52,0},{40,0}},
      color={0,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics),
Documentation(info="<html>
Example that compute the heat and mass transfer on a gypsum board with phi = 0.8, T = 5.0°C outside and phi = 0.5, T = 30°C inside in order to validate the heat and mass conduction model by comparing the results with the BuildingPysicsLibrary model (validated by the WUFI simulator). 
</html>", revisions="<html>
<ul>
<li>
Jun 19 2014, by Antoine Hudzik:<br/>
First implementation.
</li>
</ul>
</html>"));

end BuildingPhysicsLibrary;
