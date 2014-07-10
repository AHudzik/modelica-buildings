within Buildings.HeatTransfer.Conduction.BaseClasses;
model test

  SingleLayerHM lay(
    phi_ini=0.2,
    w_ini=0.003,
    activatesuction=false,
    redeclare Buildings.HeatTransfer.Data.Solids.ConcreteHM material(
      switch_lamb=1,
      switch_dw=2,
      x=0.2,
      switch_w=1),
    switch_lamb=2,
    switch_dw=2,
    T_ini=293.15,
    switch_w=3)
    annotation (Placement(transformation(extent={{-18,-24},{24,24}})));

  Sources.FixedTemperature fixedTemperature(T=288.15)
    annotation (Placement(transformation(extent={{-72,36},{-52,56}})));
  Sources.FixedTemperature fixedTemperature1(T=296.15)
    annotation (Placement(transformation(extent={{70,20},{50,40}})));
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
      points={{-18,0},{-40,0}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(lay.heatMassPort_b, heatMassPort_b1) annotation (Line(
      points={{24,0},{36,0}},
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
      points={{50,30},{44,30},{44,0},{36,0}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end test;
