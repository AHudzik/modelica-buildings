within Buildings.HeatTransfer.Conduction;
package Examples
  extends Modelica.Icons.ExamplesPackage;
  model Test
    extends Modelica.Icons.Example;

    Sources.FixedHumidity fixedHumidity(X=0.005)
      annotation (Placement(transformation(extent={{-78,-40},{-58,-20}})));

    Sources.FixedTemperature fixedTemperature1(T=293.15)
      annotation (Placement(transformation(extent={{80,20},{60,40}})));
    Sources.FixedHumidity fixedHumidity1(X=0.009)
      annotation (Placement(transformation(extent={{80,-40},{60,-20}})));
    Sources.PrescribedTemperature prescribedTemperature
      annotation (Placement(transformation(extent={{-78,14},{-58,34}})));
    Interfaces.HeatMassPort_a heatMassPort_a1
      annotation (Placement(transformation(extent={{-46,-6},{-34,6}})));
    Interfaces.HeatMassPort_b heatMassPort_b1
      annotation (Placement(transformation(extent={{34,-6},{46,6}})));
    Modelica.Blocks.Sources.Constant const(k=293.15)
      annotation (Placement(transformation(extent={{-114,14},{-94,34}})));
    SingleLayerHM singleLayerHM2_1(redeclare
        Buildings.HeatTransfer.Data.Solids.ConcreteHM material(x=0.2))
      annotation (Placement(transformation(extent={{-28,-26},{28,26}})));
  equation
    connect(prescribedTemperature.port, heatMassPort_a1.heatPort) annotation (
        Line(
        points={{-58,24},{-48,24},{-48,0},{-40,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
        points={{-58,-30},{-48,-30},{-48,0},{-40,0}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(fixedTemperature1.port, heatMassPort_b1.heatPort) annotation (Line(
        points={{60,30},{52,30},{52,0},{40,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedHumidity1.massPort, heatMassPort_b1.massPort) annotation (Line(
        points={{60,-30},{52,-30},{52,0},{40,0}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(const.y, prescribedTemperature.T) annotation (Line(
        points={{-93,24},{-80,24}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(singleLayerHM2_1.heatMassPort_a, heatMassPort_a1) annotation (Line(
        points={{-28,0},{-40,0}},
        color={0,0,0},
        pattern=LinePattern.None,
        smooth=Smooth.None));
    connect(singleLayerHM2_1.heatMassPort_b, heatMassPort_b1) annotation (Line(
        points={{28,0},{40,0}},
        color={127,0,127},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
          graphics));
  end Test;

  model TestMultiLayerHM

    Sources.PrescribedTemperature prescribedTemperature
      annotation (Placement(transformation(extent={{-72,40},{-52,60}})));
    Sources.FixedTemperature fixedTemperature(T=293.15)
      annotation (Placement(transformation(extent={{84,22},{64,42}})));
    Sources.FixedHumidity fixedHumidity(X=0.005)
      annotation (Placement(transformation(extent={{-100,-40},{-80,-20}})));
    Sources.FixedHumidity fixedHumidity1(X=0.009)
      annotation (Placement(transformation(extent={{90,-30},{70,-10}})));
    Modelica.Blocks.Sources.Constant const(k=283.15)
      annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
    Interfaces.HeatMassPort_a heatMassPort_a1
      annotation (Placement(transformation(extent={{-76,0},{-66,10}})));
    Interfaces.HeatMassPort_b heatMassPort_b1
      annotation (Placement(transformation(extent={{38,-4},{46,4}})));
    MultiLayerHM multiLayerHMbis(
      A=1,
      redeclare Buildings.HeatTransfer.Data.OpaqueConstructions.Concrete200HM
        layers,
      lay(material(Kunzel=false)))
      annotation (Placement(transformation(extent={{-36,-28},{20,28}})));
  equation
    connect(const.y, prescribedTemperature.T) annotation (Line(
        points={{-79,90},{-76,90},{-76,50},{-74,50}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
        points={{-80,-30},{-76,-30},{-76,5},{-71,5}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(prescribedTemperature.port, heatMassPort_a1.heatPort) annotation (
        Line(
        points={{-52,50},{-62,50},{-62,5},{-71,5}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedTemperature.port, heatMassPort_b1.heatPort) annotation (Line(
        points={{64,32},{64,0},{42,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedHumidity1.massPort, heatMassPort_b1.massPort) annotation (Line(
        points={{70,-20},{64,-20},{64,0},{42,0}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(heatMassPort_a1, multiLayerHMbis.heatMassPort_a) annotation (Line(
        points={{-71,5},{-53.5,5},{-53.5,0},{-36,0}},
        color={0,0,0},
        pattern=LinePattern.None,
        smooth=Smooth.None));
    connect(multiLayerHMbis.heatMassPort_b, heatMassPort_b1) annotation (Line(
        points={{20,0},{42,0}},
        color={127,0,127},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),      graphics));
  end TestMultiLayerHM;
end Examples;
