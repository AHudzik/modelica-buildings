within Buildings.HeatTransfer.Examples;
model ConductorMultiLayerHM
  extends Modelica.Icons.Example;
  Conduction.MultiLayerHM multiLayerHM(redeclare
      Buildings.HeatTransfer.Data.OpaqueConstructions.Concrete200HM layers, A=1,
    switch_w=2,
    switch_lamb=2,
    switch_dw=2,
    activatesuction=false)
    annotation (Placement(transformation(extent={{-20,-8},{16,28}})));
  Conduction.SingleLayerHM lay(redeclare
      Buildings.HeatTransfer.Data.Solids.ConcreteHM material(x=0.2, nSta=3),
    switch_w=2,
    switch_lamb=2,
    switch_dw=2,
    activatesuction=false)
    annotation (Placement(transformation(extent={{-16,-44},{14,-14}})));
  Sources.PrescribedTemperature prescribedTemperature
    annotation (Placement(transformation(extent={{-76,6},{-64,18}})));
  Sources.PrescribedTemperature prescribedTemperature1
    annotation (Placement(transformation(extent={{-76,-30},{-64,-18}})));
  Sources.FixedHumidity fixedHumidity(X=0.01)
    annotation (Placement(transformation(extent={{-76,30},{-64,42}})));
  Sources.FixedHumidity fixedHumidity1(X=0.01)
    annotation (Placement(transformation(extent={{-76,-52},{-64,-40}})));
  Interfaces.HeatMassPort_b heatMassPort_b1
    annotation (Placement(transformation(extent={{20,-32},{26,-26}})));
  Interfaces.HeatMassPort_b heatMassPort_b2
    annotation (Placement(transformation(extent={{20,6},{28,14}})));
  Sources.FixedHumidity fixedHumidity2(X=0.005)
    annotation (Placement(transformation(extent={{52,40},{40,52}})));
  Sources.FixedHumidity fixedHumidity3(X=0.005)
    annotation (Placement(transformation(extent={{52,-60},{40,-48}})));
  Modelica.Blocks.Sources.Step step(
    height=10,
    offset=293.15,
    startTime=3600)
    annotation (Placement(transformation(extent={{-100,-10},{-88,2}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T(
        displayUnit="K") = 293.15)
    annotation (Placement(transformation(extent={{80,0},{60,20}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T(
        displayUnit="K") = 293.15)
    annotation (Placement(transformation(extent={{80,-40},{60,-20}})));
  Convection.InteriorHM interiorHM(A=1, til=Buildings.HeatTransfer.Types.Tilt.Wall)
    annotation (Placement(transformation(extent={{-42,4},{-30,16}})));
  Interfaces.HeatMassPort_a solid1
    annotation (Placement(transformation(extent={{-54,8},{-50,12}})));
  Convection.InteriorHM interiorHM1(A=1, til=Buildings.HeatTransfer.Types.Tilt.Wall)
    annotation (Placement(transformation(extent={{-40,-36},{-28,-24}})));
  Interfaces.HeatMassPort_a solid2
    annotation (Placement(transformation(extent={{-52,-34},{-44,-26}})));
equation
  connect(lay.heatMassPort_b, heatMassPort_b1) annotation (Line(
      points={{14,-29},{23,-29}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(multiLayerHM.heatMassPort_b, heatMassPort_b2) annotation (Line(
      points={{16,10},{24,10}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(step.y, prescribedTemperature.T) annotation (Line(
      points={{-87.4,-4},{-82,-4},{-82,12},{-77.2,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, prescribedTemperature1.T) annotation (Line(
      points={{-87.4,-4},{-82,-4},{-82,-24},{-77.2,-24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(heatMassPort_b2.heatPort, fixedTemperature.port) annotation (Line(
      points={{24,10},{60,10}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(heatMassPort_b1.heatPort, fixedTemperature1.port) annotation (Line(
      points={{23,-29},{41.5,-29},{41.5,-30},{60,-30}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(fixedHumidity3.massPort, heatMassPort_b1.massPort) annotation (Line(
      points={{40,-54},{32,-54},{32,-29},{23,-29}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity2.massPort, heatMassPort_b2.massPort) annotation (Line(
      points={{40,46},{34,46},{34,10},{24,10}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(interiorHM.solid, solid1) annotation (Line(
      points={{-42,10.36},{-56,10.36},{-56,10},{-52,10}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(prescribedTemperature.port, solid1.heatPort) annotation (Line(
      points={{-64,12},{-60,12},{-60,10},{-52,10}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, solid1.massPort) annotation (Line(
      points={{-64,36},{-62,36},{-62,10},{-52,10}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(interiorHM.fluid, multiLayerHM.heatMassPort_a) annotation (Line(
      points={{-30,10},{-20,10}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(interiorHM1.solid, solid2) annotation (Line(
      points={{-40,-29.64},{-44,-29.64},{-44,-30},{-48,-30}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(prescribedTemperature1.port, solid2.heatPort) annotation (Line(
      points={{-64,-24},{-56,-24},{-56,-30},{-48,-30}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity1.massPort, solid2.massPort) annotation (Line(
      points={{-64,-46},{-56,-46},{-56,-30},{-48,-30}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(interiorHM1.fluid, lay.heatMassPort_a) annotation (Line(
      points={{-28,-30},{-22,-30},{-22,-29},{-16,-29}},
      color={127,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end ConductorMultiLayerHM;
