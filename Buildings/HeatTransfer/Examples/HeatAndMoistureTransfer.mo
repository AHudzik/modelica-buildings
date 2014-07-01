within Buildings.HeatTransfer.Examples;
model HeatAndMoistureTransfer
  extends Modelica.Icons.Example;
  Convection.InteriorHM interiorHM(     til=0, A=1)
    annotation (Placement(transformation(extent={{36,-8},{52,8}})));
  Convection.ExteriorHM exteriorHM(
    conMod=Buildings.HeatTransfer.Types.ExteriorConvection.Fixed,
    A=1,
    til=Buildings.HeatTransfer.Types.Tilt.Wall,
    azi=Buildings.HeatTransfer.Types.Azimuth.N)
    annotation (Placement(transformation(extent={{-42,-8},{-26,8}})));
  Modelica.Blocks.Sources.Constant vWin(k=2) "Wind speed"
    annotation (Placement(transformation(extent={{-94,58},{-78,74}})));
  Modelica.Blocks.Sources.Ramp direction1(
                                         duration=3600, height=2*3.14159)
    "Wind direction (0=from north)"
    annotation (Placement(transformation(extent={{-94,30},{-78,46}})));
  Sources.PrescribedTemperature prescribedTemperature
    annotation (Placement(transformation(extent={{-70,-20},{-56,-6}})));
  Sources.PrescribedHumidity prescribedHumidity
    annotation (Placement(transformation(extent={{-70,-40},{-56,-26}})));
  Sources.FixedTemperature fixedTemperature(T=294.15)
    annotation (Placement(transformation(extent={{100,20},{86,34}})));
  Sources.FixedHumidity fixedHumidity(X=0.0083)
    annotation (Placement(transformation(extent={{100,-34},{86,-20}})));
  Modelica.Blocks.Sources.Constant const(k=0.005)
    annotation (Placement(transformation(extent={{-100,-40},{-86,-26}})));
  Modelica.Blocks.Sources.Constant const1(k=283.15)
    annotation (Placement(transformation(extent={{-100,-18},{-86,-4}})));
  Interfaces.HeatMassPort_a solid1
    annotation (Placement(transformation(extent={{-50,-26},{-42,-18}})));
  Interfaces.HeatMassPort_b fluid1
    annotation (Placement(transformation(extent={{62,-4},{68,2}})));
  Conduction.SingleLayerHM lay(redeclare
      Buildings.HeatTransfer.Data.Solids.ConcreteHM material)
    annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
equation
  connect(fixedTemperature.port, fluid1.heatPort) annotation (Line(
      points={{86,27},{84,27},{84,-1},{65,-1}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, fluid1.massPort) annotation (Line(
      points={{86,-27},{86,-1},{65,-1}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(vWin.y, exteriorHM.v) annotation (Line(
      points={{-77.2,66},{-60,66},{-60,8},{-43.6,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(direction1.y, exteriorHM.dir) annotation (Line(
      points={{-77.2,38},{-60,38},{-60,4},{-43.6,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(exteriorHM.solid, solid1) annotation (Line(
      points={{-42,0.48},{-46,0.48},{-46,-22}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(prescribedTemperature.port, solid1.heatPort) annotation (Line(
      points={{-56,-13},{-52,-13},{-52,-22},{-46,-22}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(prescribedHumidity.massPort, solid1.massPort) annotation (Line(
      points={{-56,-32.86},{-52,-32.86},{-52,-22},{-46,-22}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(const.y, prescribedHumidity.X) annotation (Line(
      points={{-85.3,-33},{-77.65,-33},{-77.65,-33},{-71.4,-33}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const1.y, prescribedTemperature.T) annotation (Line(
      points={{-85.3,-11},{-77.65,-11},{-77.65,-13},{-71.4,-13}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(interiorHM.fluid, fluid1) annotation (Line(
      points={{52,0},{58,0},{58,-1},{65,-1}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(exteriorHM.fluid, lay.heatMassPort_a) annotation (Line(
      points={{-26,0},{-8,0}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(lay.heatMassPort_b, interiorHM.solid) annotation (Line(
      points={{12,0},{24,0},{24,0.48},{36,0.48}},
      color={127,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(coordinateSystem(extent={{-100,
            -100},{100,100}}, preserveAspectRatio=false), graphics));
end HeatAndMoistureTransfer;
