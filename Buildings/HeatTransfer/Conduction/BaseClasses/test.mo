within Buildings.HeatTransfer.Conduction.BaseClasses;
model test
  SingleLayerHM lay(redeclare Buildings.HeatTransfer.Data.Solids.ConcreteHM
      material(x=0.3, Kunzel=false))
    annotation (Placement(transformation(extent={{-20,-24},{22,24}})));
  Sources.FixedTemperature fixedTemperature(T=293.15)
    annotation (Placement(transformation(extent={{-72,36},{-52,56}})));
  Sources.FixedTemperature fixedTemperature1(T=298.15)
    annotation (Placement(transformation(extent={{70,40},{50,60}})));
  Sources.FixedHumidity fixedHumidity(X=0.005)
    annotation (Placement(transformation(extent={{-76,-26},{-56,-6}})));
  Sources.FixedHumidity fixedHumidity1(X=0.008)
    annotation (Placement(transformation(extent={{64,-32},{44,-12}})));
  Interfaces.HeatMassPort_a heatMassPort_a1
    annotation (Placement(transformation(extent={{-46,-6},{-34,6}})));
  Interfaces.HeatMassPort_b heatMassPort_b1
    annotation (Placement(transformation(extent={{32,-4},{40,4}})));
equation
  connect(lay.heatMassPort_a, heatMassPort_a1) annotation (Line(
      points={{-20,0},{-40,0}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(lay.heatMassPort_b, heatMassPort_b1) annotation (Line(
      points={{22,0},{36,0}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
      points={{-56,-16},{-48,-16},{-48,0},{-40,0}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature.port, heatMassPort_a1.heatPort) annotation (Line(
      points={{-52,46},{-46,46},{-46,0},{-40,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity1.massPort, heatMassPort_b1.massPort) annotation (Line(
      points={{44,-22},{40,-22},{40,0},{36,0}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature1.port, heatMassPort_b1.heatPort) annotation (Line(
      points={{50,50},{44,50},{44,0},{36,0}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end test;
