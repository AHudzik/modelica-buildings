within Buildings.HeatTransfer.Convection.Examples;
model InteriorHM "Test model for convective heat and mass transfer coefficient"
  import Buildings;
  extends Modelica.Icons.Example;
  Buildings.HeatTransfer.Convection.InteriorHM conCon(
    A=1,
    til=Buildings.HeatTransfer.Types.Tilt.Wall,
    conMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-50,10})));
  Buildings.HeatTransfer.Convection.InteriorHM conVer(
    A=1,
    conMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature,
    til=Buildings.HeatTransfer.Types.Tilt.Wall) annotation (Placement(
        transformation(
        extent={{10,10},{-10,-10}},
        rotation=-90,
        origin={-10,10})));
  Buildings.HeatTransfer.Convection.InteriorHM conHorFluTop(
    A=1,
    til=Buildings.HeatTransfer.Types.Tilt.Floor,
    conMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature)
    "Convection model with fluid on top" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={30,10})));
  Buildings.HeatTransfer.Convection.InteriorHM conHorSolTop(
    A=1,
    til=Buildings.HeatTransfer.Types.Tilt.Ceiling,
    conMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature)
    "Convection model with solid on top" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={70,10})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_a solid1
    annotation (Placement(transformation(extent={{-54,-26},{-44,-16}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_a solid2
    annotation (Placement(transformation(extent={{-16,-26},{-6,-16}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_a solid3
    annotation (Placement(transformation(extent={{24,-26},{34,-16}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_a solid4
    annotation (Placement(transformation(extent={{64,-26},{74,-16}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_b fluid1
    annotation (Placement(transformation(extent={{-54,36},{-44,46}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_b fluid2
    annotation (Placement(transformation(extent={{-16,36},{-6,46}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_b fluid3
    annotation (Placement(transformation(extent={{24,36},{34,46}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_b fluid4
    annotation (Placement(transformation(extent={{64,36},{74,46}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature TA1
    annotation (Placement(transformation(extent={{-72,60},{-60,72}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature TA2
    annotation (Placement(transformation(extent={{-32,60},{-20,72}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature TA3
    annotation (Placement(transformation(extent={{8,60},{20,72}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature TA4
    annotation (Placement(transformation(extent={{48,60},{60,72}})));
  Buildings.HeatTransfer.Sources.PrescribedHumidity HA1
    annotation (Placement(transformation(extent={{-72,80},{-60,92}})));
  Buildings.HeatTransfer.Sources.PrescribedHumidity HA2
    annotation (Placement(transformation(extent={{-32,80},{-20,92}})));
  Buildings.HeatTransfer.Sources.PrescribedHumidity HA3
    annotation (Placement(transformation(extent={{8,80},{20,92}})));
  Buildings.HeatTransfer.Sources.PrescribedHumidity HA4
    annotation (Placement(transformation(extent={{48,80},{60,92}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TB1(T=293.15)
    annotation (Placement(transformation(extent={{-72,-52},{-60,-40}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TB2(T=293.15)
    annotation (Placement(transformation(extent={{-30,-52},{-18,-40}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TB3(T=293.15)
    annotation (Placement(transformation(extent={{8,-52},{20,-40}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TB4(T=293.15)
    annotation (Placement(transformation(extent={{48,-52},{60,-40}})));
  Buildings.HeatTransfer.Sources.FixedHumidity HB1(X=0.0115)
    annotation (Placement(transformation(extent={{-72,-72},{-60,-60}})));
  Buildings.HeatTransfer.Sources.FixedHumidity HB2(X=0.0115)
    annotation (Placement(transformation(extent={{-32,-72},{-20,-60}})));
  Buildings.HeatTransfer.Sources.FixedHumidity HB3(X=0.0115)
    annotation (Placement(transformation(extent={{8,-72},{20,-60}})));
  Buildings.HeatTransfer.Sources.FixedHumidity HB4(X=0.0115)
    annotation (Placement(transformation(extent={{48,-72},{60,-60}})));
  Modelica.Blocks.Sources.Ramp TempStep(
    height=10,
    duration=1,
    offset=293.15 - 5)
    annotation (Placement(transformation(extent={{-118,44},{-98,64}})));
  Modelica.Blocks.Sources.Ramp HumidityStep(
    height=0.01,
    duration=1,
    offset=0.005)
    annotation (Placement(transformation(extent={{-118,90},{-98,110}})));
equation
  connect(conCon.solid, solid1) annotation (Line(
      points={{-50.6,0},{-49,0},{-49,-21}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(conVer.solid, solid2) annotation (Line(
      points={{-10.6,0},{-11,0},{-11,-21}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(conHorFluTop.solid, solid3) annotation (Line(
      points={{29.4,0},{29,0},{29,-21}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(conHorSolTop.solid, solid4) annotation (Line(
      points={{69.4,0},{70,0},{70,-21},{69,-21}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(conCon.fluid, fluid1) annotation (Line(
      points={{-50,20},{-50,41},{-49,41}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(conVer.fluid, fluid2) annotation (Line(
      points={{-10,20},{-10,41},{-11,41}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(conHorFluTop.fluid, fluid3) annotation (Line(
      points={{30,20},{30,41},{29,41}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(conHorSolTop.fluid, fluid4) annotation (Line(
      points={{70,20},{70,41},{69,41}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(TempStep.y, TA1.T) annotation (Line(
      points={{-97,54},{-80,54},{-80,66},{-73.2,66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TempStep.y, TA2.T) annotation (Line(
      points={{-97,54},{-40,54},{-40,66},{-33.2,66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TempStep.y, TA3.T) annotation (Line(
      points={{-97,54},{0,54},{0,66},{6.8,66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TempStep.y, TA4.T) annotation (Line(
      points={{-97,54},{40,54},{40,66},{46.8,66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(HumidityStep.y, HA1.X) annotation (Line(
      points={{-97,100},{-84,100},{-84,86},{-73.2,86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(HumidityStep.y, HA2.X) annotation (Line(
      points={{-97,100},{-40,100},{-40,86},{-33.2,86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(HumidityStep.y, HA3.X) annotation (Line(
      points={{-97,100},{0,100},{0,86},{6.8,86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(HumidityStep.y, HA4.X) annotation (Line(
      points={{-97,100},{40,100},{40,86},{46.8,86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TA1.port, fluid1.heatPort) annotation (Line(
      points={{-60,66},{-58,66},{-58,41},{-49,41}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(HA1.massPort, fluid1.massPort) annotation (Line(
      points={{-60,86.12},{-50,86.12},{-50,41},{-49,41}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(TA2.port, fluid2.heatPort) annotation (Line(
      points={{-20,66},{-18,66},{-18,41},{-11,41}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(TA3.port, fluid3.heatPort) annotation (Line(
      points={{20,66},{22,66},{22,41},{29,41}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(TA4.port, fluid4.heatPort) annotation (Line(
      points={{60,66},{62,66},{62,41},{69,41}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(HA2.massPort, fluid2.massPort) annotation (Line(
      points={{-20,86.12},{-12,86.12},{-12,41},{-11,41}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(HA3.massPort, fluid3.massPort) annotation (Line(
      points={{20,86.12},{28,86.12},{28,41},{29,41}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(HA4.massPort, fluid4.massPort) annotation (Line(
      points={{60,86.12},{68,86.12},{68,41},{69,41}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(TB1.port, solid1.heatPort) annotation (Line(
      points={{-60,-46},{-58,-46},{-58,-21},{-49,-21}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(TB2.port, solid2.heatPort) annotation (Line(
      points={{-18,-46},{-18,-21},{-11,-21}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(TB3.port, solid3.heatPort) annotation (Line(
      points={{20,-46},{22,-46},{22,-21},{29,-21}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(TB4.port, solid4.heatPort) annotation (Line(
      points={{60,-46},{62,-46},{62,-21},{69,-21}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(HB1.massPort, solid1.massPort) annotation (Line(
      points={{-60,-66},{-50,-66},{-50,-21},{-49,-21}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(HB2.massPort, solid2.massPort) annotation (Line(
      points={{-20,-66},{-12,-66},{-12,-21},{-11,-21}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(HB3.massPort, solid3.massPort) annotation (Line(
      points={{20,-66},{28,-66},{28,-21},{29,-21}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(HB4.massPort, solid4.massPort) annotation (Line(
      points={{60,-66},{68,-66},{68,-21},{69,-21}},
      color={0,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,
            -120},{120,120}}), graphics), Icon(coordinateSystem(extent={{-120,
            -120},{120,120}})));
end InteriorHM;
