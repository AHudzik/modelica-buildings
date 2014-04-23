within Buildings.HeatTransfer.Conduction;
package Examples
  extends Modelica.Icons.ExamplesPackage;
  model Test
    extends Modelica.Icons.Example;

    Sources.FixedHumidity fixedHumidity(X=0.009)
      annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

    Sources.FixedTemperature fixedTemperature1(T=308.15)
      annotation (Placement(transformation(extent={{80,20},{60,40}})));
    Sources.FixedHumidity fixedHumidity1(X=0.0164)
      annotation (Placement(transformation(extent={{80,-40},{60,-20}})));
    Sources.PrescribedTemperature prescribedTemperature
      annotation (Placement(transformation(extent={{-78,14},{-58,34}})));
    Interfaces.HeatMassPort_a heatMassPort_a1
      annotation (Placement(transformation(extent={{-46,-6},{-34,6}})));
    Interfaces.HeatMassPort_b heatMassPort_b1
      annotation (Placement(transformation(extent={{34,-6},{46,6}})));
    Modelica.Blocks.Sources.Constant const(k=298.15)
      annotation (Placement(transformation(extent={{-114,14},{-94,34}})));
    SingleLayerHM singleLayerHM2_1
      annotation (Placement(transformation(extent={{-28,-26},{28,26}})));
  equation
    connect(prescribedTemperature.port, heatMassPort_a1.heatPort) annotation (
        Line(
        points={{-58,24},{-48,24},{-48,0},{-40,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
        points={{-60,-30},{-48,-30},{-48,0},{-40,0}},
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
end Examples;
