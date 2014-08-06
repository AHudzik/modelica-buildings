within Buildings.HeatTransfer.Examples;
model ConductorSingleLayerHM
  extends Modelica.Icons.Example;
  Modelica.Blocks.Sources.Step step(
    height=10,
    offset=293.15,
    startTime=3600)
    annotation (Placement(transformation(extent={{-100,-26},{-88,-14}})));
  Interfaces.HeatMassPort_b heatMassPort_b2
    annotation (Placement(transformation(extent={{16,-44},{24,-36}})));
  Sources.FixedHumidity fixedHumidity2(X=0.005)
    annotation (Placement(transformation(extent={{76,-20},{64,-8}})));
  Sources.FixedHumidity fixedHumidity3(X=0.005)
    annotation (Placement(transformation(extent={{64,-60},{52,-48}})));
  Sources.FixedTemperature fixedTemperature(T=293.15)
    annotation (Placement(transformation(extent={{92,-6},{80,6}})));
  Sources.FixedTemperature fixedTemperature1(T=293.15)
    annotation (Placement(transformation(extent={{72,-46},{60,-34}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-76,-6},{-64,6}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature1
    annotation (Placement(transformation(extent={{-76,-46},{-64,-34}})));
  Sources.FixedHumidity fixedHumidity(X=0.01)
    annotation (Placement(transformation(extent={{-66,20},{-56,30}})));
  Sources.FixedHumidity fixedHumidity1(X=0.01)
    annotation (Placement(transformation(extent={{-66,-70},{-56,-60}})));
  Interfaces.HeatMassPort_b heatMassPort_b1
    annotation (Placement(transformation(extent={{36,-4},{44,4}})));
  Conduction.SingleLayerHM lay1(redeclare
      Buildings.HeatTransfer.Data.Solids.ConcreteHM material(x=0.1, Kunzel=
          false))
    annotation (Placement(transformation(extent={{-14,-10},{6,10}})));
  Conduction.SingleLayerHM lay2(redeclare
      Buildings.HeatTransfer.Data.Solids.ConcreteHM material(x=0.1, Kunzel=
          false))
    annotation (Placement(transformation(extent={{10,-10},{30,10}})));
  Conduction.MultiLayerHM multiLayerHMbis(redeclare
      Buildings.HeatTransfer.Data.OpaqueConstructions.Concrete200HM layers, A=1)
    annotation (Placement(transformation(extent={{-8,-50},{12,-30}})));
  Convection.InteriorHM interiorHM(A=1, til=0)
    annotation (Placement(transformation(extent={{-32,-6},{-20,6}})));
  Convection.InteriorHM interiorHM1(A=1, til=0)
    annotation (Placement(transformation(extent={{-32,-46},{-20,-34}})));
  Interfaces.HeatMassPort_a solid1
    annotation (Placement(transformation(extent={{-48,-4},{-40,4}})));
  Interfaces.HeatMassPort_a solid2
    annotation (Placement(transformation(extent={{-46,-44},{-38,-36}})));
equation
  connect(fixedTemperature1.port, heatMassPort_b2.heatPort) annotation (Line(
      points={{60,-40},{20,-40}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(step.y, prescribedTemperature.T) annotation (Line(
      points={{-87.4,-20},{-82,-20},{-82,0},{-77.2,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, prescribedTemperature1.T) annotation (Line(
      points={{-87.4,-20},{-82,-20},{-82,-40},{-77.2,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(fixedHumidity3.massPort, heatMassPort_b2.massPort) annotation (Line(
      points={{52,-54},{36,-54},{36,-40},{20,-40}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity2.massPort, heatMassPort_b1.massPort) annotation (Line(
      points={{64,-14},{52,-14},{52,0},{40,0}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature.port, heatMassPort_b1.heatPort) annotation (Line(
      points={{80,0},{40,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(lay1.heatMassPort_b, lay2.heatMassPort_a) annotation (Line(
      points={{6,0},{10,0}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(lay2.heatMassPort_b, heatMassPort_b1) annotation (Line(
      points={{30,0},{40,0}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(multiLayerHMbis.heatMassPort_b, heatMassPort_b2) annotation (Line(
      points={{12,-40},{20,-40}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(interiorHM.solid, solid1) annotation (Line(
      points={{-32,0.36},{-38,0.36},{-38,0},{-44,0}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(interiorHM1.solid, solid2) annotation (Line(
      points={{-32,-39.64},{-38,-39.64},{-38,-40},{-42,-40}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(prescribedTemperature.port, solid1.heatPort) annotation (Line(
      points={{-64,0},{-44,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, solid1.massPort) annotation (Line(
      points={{-56,25},{-50,25},{-50,0},{-44,0}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(prescribedTemperature1.port, solid2.heatPort) annotation (Line(
      points={{-64,-40},{-42,-40}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity1.massPort, solid2.massPort) annotation (Line(
      points={{-56,-65},{-50,-65},{-50,-40},{-42,-40}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(interiorHM.fluid, lay1.heatMassPort_a) annotation (Line(
      points={{-20,0},{-14,0}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(interiorHM1.fluid, multiLayerHMbis.heatMassPort_a) annotation (Line(
      points={{-20,-40},{-8,-40}},
      color={127,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end ConductorSingleLayerHM;